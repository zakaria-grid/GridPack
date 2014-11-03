/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * PTI23parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef PARSER_HPP_
#define PARSER_HPP_

#define OLD_MAP

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "gridpack/utilities/exception.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class _network>
class Parser
{
    public:
  /// Constructor 
  /**
   * 
   * @param network network object that will be filled with contents
   * of network configuration file (must be child of network::BaseNetwork<>)
   */
  Parser(boost::shared_ptr<_network> network)
    : p_network(network), p_configExists(false)
  {
  }

      /**
       * Destructor
       */
      virtual ~Parser()
      {
        p_busData.clear();      // unnecessary
        p_branchData.clear();
      }

      /**
       * Parse network configuration file and create network
       * @param fileName name of network file
       */
      void parse(const std::string &fileName)
      {
          p_timer = gridpack::utility::CoarseTimer::instance();
          p_timer->configTimer(false);
          int t_total = p_timer->createCategory("Parser:Total Elapsed Time");
          p_timer->start(t_total);

          gridpack::parser::GOSSParser           parser;

          try {
              parser.parse(fileName.c_str());
          } catch (std::exception & e) {
              e.what();
          } catch (boost::exception & e) {
              std::cout << "General exception\n\t";
              std::cout << std::endl;
          }

          p_case_id      = parser.getCaseId();
          p_case_sbase   = parser.getCaseSbase();
          parser.loadBuses(p_busData,p_busMap);
          parser.loadBranches(p_branchData);

          createNetwork();

          p_timer->stop(t_total);
          p_timer->configTimer(true);

      }

    protected:

      void createNetwork(void)
      {
        int t_create = p_timer->createCategory("Parser:createNetwork");
        p_timer->start(t_create);
        int me(p_network->communicator().rank());
        int nprocs(p_network->communicator().size());
        int i;
        // Exchange information on number of buses and branches on each
        // processor
        int sbus[nprocs], sbranch[nprocs];
        int nbus[nprocs], nbranch[nprocs];
        for (i=0; i<nprocs; i++) {
          sbus[i] = 0;
          sbranch[i] = 0;
        }
        sbus[me] = p_busData.size();
        sbranch[me] = p_branchData.size();
        MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
        int ierr;
        ierr = MPI_Allreduce(sbus,nbus,nprocs,MPI_INT,MPI_SUM,comm);
        ierr = MPI_Allreduce(sbranch,nbranch,nprocs,MPI_INT,MPI_SUM,comm);
        // Transmit CASE_ID and CASE_SBASE to all processors
        int isval, irval;
        if (me == 0) {
          isval = p_case_id;
        } else {
          isval = 0;
        }
        ierr = MPI_Allreduce(&isval,&irval,1,MPI_INT,MPI_SUM,comm);
        p_case_id = irval;
        double sval, rval;
        if (me == 0) {
          sval = p_case_sbase;
        } else {
          sval = 0.0;
        }
        ierr = MPI_Allreduce(&sval,&rval,1,MPI_DOUBLE,MPI_SUM,comm);
        p_case_sbase = rval;
        // evaluate offsets for buses and branches
        int offset_bus[nprocs], offset_branch[nprocs];
        offset_bus[0] = 0;
        offset_branch[0] = 0;
        for (i=1; i<nprocs; i++) {
          offset_bus[i] = offset_bus[i-1]+nbus[i-1];
          offset_branch[i] = offset_branch[i-1]+nbranch[i-1];
        }

        int numBus = p_busData.size();
        for (i=0; i<numBus; i++) {
          int idx;
          p_busData[i]->getValue(BUS_NUMBER,&idx);
          p_network->addBus(idx);
          p_network->setGlobalBusIndex(i,i+offset_bus[me]);
          *(p_network->getBusData(i)) = *(p_busData[i]);
          p_network->getBusData(i)->addValue(CASE_ID,p_case_id);
          p_network->getBusData(i)->addValue(CASE_SBASE,p_case_sbase);
        }
        int numBranch = p_branchData.size();
        for (i=0; i<numBranch; i++) {
          int idx1, idx2;
          p_branchData[i]->getValue(BRANCH_FROMBUS,&idx1);
          p_branchData[i]->getValue(BRANCH_TOBUS,&idx2);
          p_network->addBranch(idx1, idx2);
          p_network->setGlobalBranchIndex(i,i+offset_branch[me]);
          int g_idx1, g_idx2;
          std::map<int, int>::iterator it;

          it = p_busMap.find(idx1);
          g_idx1 = it->second;
          it = p_busMap.find(idx2);
          g_idx2 = it->second;
          *(p_network->getBranchData(i)) = *(p_branchData[i]);
          p_network->getBranchData(i)->addValue(CASE_ID,p_case_id);
          p_network->getBranchData(i)->addValue(CASE_SBASE,p_case_sbase);
        }
        p_configExists = true;
        printf("Number of buses: %d\n",p_network->numBuses());
        p_busData.clear();
        p_branchData.clear();
        p_timer->stop(t_create);
      }

      /**
       * This routine opens up a .dyr file with parameters for dynamic
       * simulation. It assumes that a .raw file has already been parsed
       */
      void getDS(const std::string & fileName)
      {

        if (!p_configExists) return;
        int t_ds = p_timer->createCategory("Parser:getDS");
        p_timer->start(t_ds);
        int me(p_network->communicator().rank());

        if (me == 0) {
          std::ifstream            input;
          input.open(fileName.c_str());
          if (!input.is_open()) {
            p_timer->stop(t_ds);
            return;
          }
          find_ds_par(input);
          input.close();
        }
        p_timer->stop(t_ds);
      }

      // Clean up 2 character tags so that single quotes are removed and single
      // character tags are right-justified. These tags can be delimited by a
      // pair of single quotes, a pair of double quotes, or no quotes

      // Distribute data uniformly on processors
      void brdcst_data(void)
      {
        int t_brdcst = p_timer->createCategory("Parser:brdcst_data");
        int t_serial = p_timer->createCategory("Parser:data packing and unpacking");
        MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
        int me(p_network->communicator().rank());
        int nprocs(p_network->communicator().size());
        if (nprocs == 1) return;
        p_timer->start(t_brdcst);

        // find number of buses and branches and broadcast this information to
        // all processors
        int sbus, sbranch;
        if (me == 0) {
          sbus = p_busData.size();
          sbranch = p_branchData.size();
        } else {
          sbus = 0;
          sbranch = 0;
        }
        int ierr, nbus, nbranch;
        ierr = MPI_Allreduce(&sbus, &nbus, 1, MPI_INT, MPI_SUM, comm);
        ierr = MPI_Allreduce(&sbranch, &nbranch, 1, MPI_INT, MPI_SUM, comm);
        double rprocs = static_cast<double>(nprocs);
        double rme = static_cast<double>(me);
        int n, i;
        std::vector<gridpack::component::DataCollection>recvV;
        // distribute buses
        if (me == 0) {
          for (n=0; n<nprocs; n++) {
            double rn = static_cast<double>(n);
            int istart = static_cast<int>(static_cast<double>(nbus)*rn/rprocs);
            int iend = static_cast<int>(static_cast<double>(nbus)*(rn+1.0)/rprocs);
            if (n != 0) {
              p_timer->start(t_serial);
              std::vector<gridpack::component::DataCollection> sendV;
              for (i=istart; i<iend; i++) {
                sendV.push_back(*(p_busData[i]));
              }
              p_timer->stop(t_serial);
              static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
            } else {
              p_timer->start(t_serial);
              for (i=istart; i<iend; i++) {
                recvV.push_back(*(p_busData[i]));
              }
              p_timer->stop(t_serial);
            }
          }
        } else {
          int istart = static_cast<int>(static_cast<double>(nbus)*rme/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbus)*(rme+1.0)/rprocs)-1;
          static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
        }
        int nsize = recvV.size();
        p_busData.clear();
        p_timer->start(t_serial);
        for (i=0; i<nsize; i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data(new
              gridpack::component::DataCollection);
          *data = recvV[i];
          p_busData.push_back(data);
        }
        p_timer->stop(t_serial);
        recvV.clear();
        // distribute branches
        if (me == 0) {
          for (n=0; n<nprocs; n++) {
            double rn = static_cast<double>(n);
            int istart = static_cast<int>(static_cast<double>(nbranch)*rn/rprocs);
            int iend = static_cast<int>(static_cast<double>(nbranch)*(rn+1.0)/rprocs);
            if (n != 0) {
              p_timer->start(t_serial);
              std::vector<gridpack::component::DataCollection> sendV;
              for (i=istart; i<iend; i++) {
                sendV.push_back(*(p_branchData[i]));
              }
              p_timer->stop(t_serial);
              static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
            } else {
              p_timer->start(t_serial);
              for (i=istart; i<iend; i++) {
                recvV.push_back(*(p_branchData[i]));
              }
              p_timer->stop(t_serial);
            }
          }
        } else {
          int istart = static_cast<int>(static_cast<double>(nbranch)*rme/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbranch)*(rme+1.0)/rprocs)-1;
          static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
        }
        nsize = recvV.size();
        p_branchData.clear();
        p_timer->start(t_serial);
        for (i=0; i<nsize; i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data(new
              gridpack::component::DataCollection);
          *data = recvV[i];
          p_branchData.push_back(data);
        }
        p_timer->stop(t_serial);
        p_timer->stop(t_brdcst);
      }

    private:
     /*
       * The case_data is the collection of all data points in the case file.
       * Each collection in the case data contains the data associated with a given
       * type. For example, the case is the collection of data describing the
       * current case and the bus data is the collection of data associated with
       * each bus. The type data may consist of zero or more instances of the
       * given type. For example, the bus data may contain several instances of
       * a bus. These type instances are composed of a set of key value pairs.
       * Each column as an associated key and each row is an instance of a given
       * type. When the parser is reading data for a type, the value found in each
       * column associated with the key for that column in a field_data structure.
       *
       * Within the PTI file there are the following group of data sets in order:
       *     case
       *     bus
       *     generator
       *     branch
       *     transformer
       *     dc_line
       *     shunt
       *     impedence corr
       *     multi-terminal
       *     multi-section
       *     zone
       *     inter-area
       *     owner
       *     device driver
       *
       * These data sets are stored in the case data as a collection of
       * data set and each data set is a
       */
      boost::shared_ptr<_network> p_network;

      bool p_configExists;

      // Vector of bus data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
      // Vector of branch data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
      // Map of PTI indices to index in p_busData
      std::map<int,int> p_busMap;

      // Global variables that apply to whole network
      int p_case_id;
      double p_case_sbase;
      gridpack::utility::CoarseTimer *p_timer;

};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PARSER_HPP_ */
