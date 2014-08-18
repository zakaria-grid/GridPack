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

#ifndef PTI23_PARSER_HPP_
#define PTI23_PARSER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>

#define OLD_MAP

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif
#include <string>
#include <map>
#include <vector>
#include <exception>
#include <iostream>
#include <cstring>
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

// ----------------------------------------------------------------
// XML CODE
// ----------------------------------------------------------------
enum XML_TYPE {BOOLEAN, INTEGER, DOUBLE, CHARACTER, STRING};

template <class _network>
class PTI23_parser
{
public:
    /// Constructor
    /**
     *
     * @param network network object that will be filled with contents
     * of network configuration file (must be child of network::BaseNetwork<>)
     */
    explicit PTI23_parser(boost::shared_ptr<_network> network) : p_network(network){ }

    /**
     * Destructor
     */
    virtual ~PTI23_parser()
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
        getCase(fileName);
        //brdcst_data();
        createNetwork();
        p_timer->stop(t_total);
        p_timer->configTimer(true);
    }

    /*
     * A case is the collection of all data associated with a PTI23 file.
     * Each case is a a vector of data_set objects the contain all the data
     * associated with a partition of the PTI file. For example, the bus
     * data in the file constitutes a data_set. Each data_set is a vector of
     * gridpack::component::DataCollection objects. Each of these objects
     * contain a single instance of the data associated with a data_set. For
     * example, each line of the bus partition corresponds to a single
     * DataCollection object.
     */
    void getCase(const std::string & xmlFileName)
    {

        int t_case = p_timer->createCategory("Parser:getCase");
        p_timer->start(t_case);
        p_busData.clear();
        p_branchData.clear();
        p_busMap.clear();

        int me(p_network->communicator().rank());

        if (me == 0) {
            // ----------------------------------------------------------------
            // XML CODE
            // ----------------------------------------------------------------

            boost::property_tree::ptree xmlTree;
            read_xml(xmlFileName, xmlTree);
            std::string name("");

            find_case();
            boost::property_tree::ptree schemaStyleAttr = xmlTree.get_child("application.grammars.xs:schema");

            BOOST_FOREACH( boost::property_tree::ptree::value_type const& elementAttr, schemaStyleAttr)
            {
                setTypeAssociations(elementAttr, name);
            }

            // NOTE: the GOSS file does not have case information
            boost::property_tree::ptree branchesXMLTree = xmlTree.get_child("application.GridpackPowergrid.Branches");

            int                  index = 0;
            BOOST_FOREACH( boost::property_tree::ptree::value_type const& branchAttr, branchesXMLTree)
            {
                readBranch(branchAttr, index);
                ++ index;
            }
            boost::property_tree::ptree busXMLTree = xmlTree.get_child("application.GridpackPowergrid.Buses");

            index = 0;
            BOOST_FOREACH( boost::property_tree::ptree::value_type const& busAttr, busXMLTree)
            {
                readBus(busAttr, index);
                ++ index;
            }
            // ----------------------------------------------------------------
            // END XML CODE
            // ----------------------------------------------------------------
        }
        p_timer->stop(t_case);
    }
#if 0
          find_imped_corr(input);
          find_multi_term(input);
          find_multi_section(input);
          find_zone(input);
          find_interarea(input);
          find_owner(input);
          //find_line(input);
#endif
#if 0
          // debug
          int i;
          printf("BUS data size: %d\n",p_busData.size());
          for (i=0; i<p_busData.size(); i++) {
          printf("Dumping bus: %d\n",i);
            p_busData[i]->dump();
          }
          printf("BRANCH data size: %d\n",p_branchData.size());
          for (i=0; i<p_branchData.size(); i++) {
            printf("Dumping branch: %d\n",i);
            p_branchData[i]->dump();
          }
#endif

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
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap.find(idx1);
          g_idx1 = it->second;
          it = p_busMap.find(idx2);
          g_idx2 = it->second;
          *(p_network->getBranchData(i)) = *(p_branchData[i]);
          p_network->getBranchData(i)->addValue(CASE_ID,p_case_id);
          p_network->getBranchData(i)->addValue(CASE_SBASE,p_case_sbase);
        }
#if 0
        // debug
        printf("Number of buses: %d\n",numBus);
        for (i=0; i<numBus; i++) {
          printf("Dumping bus: %d\n",i);
          p_network->getBusData(i)->dump();
        }
        printf("Number of branches: %d\n",numBranch);
        for (i=0; i<numBranch; i++) {
          printf("Dumping branch: %d\n",i);
          p_network->getBranchData(i)->dump();
        }
#endif
        p_busData.clear();
        p_branchData.clear();
        p_timer->stop(t_create);
      }


std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & getPBusData(){return p_busData;};
std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & getPBranchData(){return p_branchData;};

protected:

      void find_case(std::ifstream & input)
      {
  //      data_set                                           case_set;
        std::string                                        line;
  //      std::vector<gridpack::component::DataCollection>   case_instance;

  //      gridpack::component::DataCollection                data;

        std::getline(input, line);
        while (check_comment(line)) {
          std::getline(input, line);
        }
        std::vector<std::string>  split_line;

        boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

        // CASE_ID             "IC"                   ranged integer
        p_case_id = 1;//atoi(split_line[0].c_str());

        // CASE_SBASE          "SBASE"                float
        p_case_sbase = 1.0; //atof(split_line[1].c_str());

        /*  These do not appear in the dictionary
        // CASE_RECORD2        "RECORD2"              string
        std::getline(input, line);
        data.addValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        // CASE_RECORD3        "RECORD3"              string
        std::getline(input, line);
        data.addValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);
         */
 //       case_set.push_back(case_instance);
 //       case_data->push_back(case_set);

      }
      //  ---------------------------------------------------------------------

// ----------------------------------------------------------------
// XML CODE
//   ----------------------------------------------------------------

/*
 * Gridpack needs the data values for p_case_id and p_case_sbase,
 * but they are not currently in the XML file. These variables will
 * be set to dummy values until the problem is resolved.
 */
void find_case()
{
    //
    // CASE_ID             "IC"                   ranged integer
    p_case_id = 1;//atoi(split_line[0].c_str());

    // CASE_SBASE          "SBASE"                float
    p_case_sbase = 1.0; //atof(split_line[1].c_str());

    /*  These do not appear in the dictionary
     * // CASE_RECORD2        "RECORD2"              string
     * std::getline(input, line);
     * data.addValue(CASE_RECORD2, line.c_str());
     * case_instance.push_back(data);
     *
     * // CASE_RECORD3        "RECORD3"              string
     * std::getline(input, line);
     * data.addValue(CASE_RECORD3, line.c_str());
     * case_instance.push_back(data);
     */
    //       case_set.push_back(case_instance);
    //       case_data->push_back(case_set);

      }

void setTypeAssociations(boost::property_tree::ptree::value_type const& elementAttr, std::string & name)
{
    XML_TYPE typeEnum;
    if (elementAttr.first == "xs:element") {
        const boost::property_tree::ptree & attributes = elementAttr.second.get_child("<xmlattr>");
        BOOST_FOREACH( boost::property_tree::ptree::value_type const& attr, attributes)
        {
            std::string               type   = attr.second.data();
            std::vector<std::string>  split_type;
            boost::algorithm::split(split_type, type, boost::algorithm::is_any_of(":"), boost::token_compress_on);
            typeMap["mrid"] = STRING;
            typeMap["elementIndex"] = INTEGER;
            if (split_type[0] == "xs"){
                if (split_type[1] == "int") {
                    typeEnum = INTEGER;
                } else if (split_type[1] == "boolean") {
                     typeEnum = BOOLEAN;
                } else if (split_type[1] == "double") {
                    typeEnum = DOUBLE;
                } else if (split_type[1] == "char") {
                    typeEnum = CHARACTER;
                } else if (split_type[1] == "string") {
                    typeEnum = STRING;
                } else {
                    // TODO: exception, invalid type
                }

                typeMap[name] = typeEnum;
            } else {
                name = attr.second.data();
            }
        }
    } else if (elementAttr.first == "xs:complexType") {
        // a complex type is a compound data type in the XML file
        boost::property_tree::ptree complexTree = elementAttr.second.get_child("xs:sequence");
        BOOST_FOREACH( boost::property_tree::ptree::value_type const& sequenceAttr, complexTree)
        {
            setTypeAssociations(sequenceAttr, name);
        }
    }
}

void loadCollection(boost::shared_ptr<gridpack::component::DataCollection> & data,
        boost::property_tree::ptree::value_type & attr)
{
    XML_TYPE         type               = typeMap[attr.first];
    switch (type)
    {
        case BOOLEAN:
        {
            bool             value          = true;
            if (attr.second.data() == "false") value = false;
            data->addValue(attr.first.c_str(), value);
            break;
        }
        case INTEGER:
        {
              data->addValue(attr.first.c_str(), atoi(attr.second.data().c_str()));
            break;
        }
        case DOUBLE:
        {
              data->addValue(attr.first.c_str(), atof(attr.second.data().c_str()));
            break;
        }
        case CHARACTER:
        {
            data->addValue(attr.first.c_str(), attr.second.data()[0]);
            break;
        }
        case STRING:
        {
            data->addValue(attr.first.c_str(), attr.second.data().c_str());
            break;
        }
    }
}

void loadCollection(boost::shared_ptr<gridpack::component::DataCollection> & data,
        boost::property_tree::ptree::value_type & attr, int nElems)
{
    XML_TYPE         type               = typeMap[attr.first];
    switch (type)
    {
        case BOOLEAN:
        {
            bool             value          = true;
            if (attr.second.data() == "false") value = false;
            data->addValue(attr.first.c_str(), value, nElems);
            break;
        }
        case INTEGER:
        {
            data->addValue(attr.first.c_str(), atoi(attr.second.data().c_str()), nElems);
            break;
        }
        case DOUBLE:
        {
            data->addValue(attr.first.c_str(), atof(attr.second.data().c_str()), nElems);
            break;
        }
        case CHARACTER:
        {
            data->addValue(attr.first.c_str(), attr.second.data()[0], nElems);
            break;
        }
        case STRING:
        {
            data->addValue(attr.first.c_str(), attr.second.data().c_str(), nElems);
            break;
        }
    }
}

void readBranch(boost::property_tree::ptree::value_type const& branchesAttr, int index)
{
    boost::shared_ptr<gridpack::component::DataCollection>  data;
    int                      connectionMade = 0;

    BOOST_FOREACH( boost::property_tree::ptree::value_type branchAttr, branchesAttr.second)
    {
        if (branchAttr.first == "BRANCH_FROMBUS") {
            data->addValue(branchAttr.first.c_str(), atoi(branchAttr.second.data().c_str()));
            ++connectionMade;
        } else  if (branchAttr.first == "BRANCH_TOBUS") {
            data->addValue(branchAttr.first.c_str(), atoi(branchAttr.second.data().c_str()));
            ++connectionMade;
        } else {
            // both BRANCH_FROMBUS and BRANCH_TOBUS must be set
            if (connectionMade != 2) {
                // TODO: throw an exception, the necessary values were not provided
            }

            // setup from and to branch relationship
            std::pair<int, int> branch_pair;
            int                 branch_fromBus      = 0;
            int                 branch_toBus        = 0;
            data->getValue(BRANCH_FROMBUS,&branch_fromBus);
            data->getValue(BRANCH_TOBUS,&branch_toBus);

            if (branch_fromBus < 0) branch_fromBus = -branch_fromBus;
            if (branch_toBus < 0) branch_toBus = -branch_toBus;
            branch_pair = std::pair<int,int>(branch_fromBus, branch_toBus);

            // setup iterator used to retrieve the branch pair if it is stored in the map
        #ifdef OLD_MAP
            std::map<std::pair<int, int>, int>::iterator it;
        #else
            boost::unordered_map<std::pair<int, int>, int>::iterator it;
        #endif

            it = p_branchMap.find(branch_pair);

            // The branch map contains the following
            //      - pair consisting of from and to bus connections
            //      - index, which increments for each bus read from the file
            //
            // On the first occurrence of an edge between two buses:
            //   A vertex pair (from and to bus indices) is mapped to an index, the
            //      index is the size of the branch map.
            //      NOTE: the edge can be a->b or b->a
            //   The "branch_index" (l_idx) is set to the size of branch data vector
            //   A branch data object is created and inserted in to the branch data
            //      vector
            // If the vertex pair is found in the branch map
            //   The branch_index is set to the corresponding map index (the index of
            //      the mapping within the branch map.
            //   The value of nelems is retrieved from the branch data at branch_index
            //      nelems is the index of the branch data within a given edge
            // If the vertex pair is not in the branch map
            //   Search for switched branch (to and from) in the branch map
            //   branch_index is the index in the branch map corresponding the first
            //      occurrence of the element in the branch map
            //
            // For each branch model parameter, the values are a triple of the name of
            //  the parameter, the value of the parameter and the number of elements in
            //  the parameter.
            bool                     switched       = false;

            int                      branch_index   = 0;
            int                      nElems         = 0;
            if (it != p_branchMap.end()) {
                branch_index = it->second;
                p_branchData[branch_index]->getValue(BRANCH_NUM_ELEMENTS, &nElems);
            } else {
                std::pair<int, int> new_branch_pair;
                new_branch_pair = std::pair<int,int>(branch_toBus, branch_fromBus);
                it = p_branchMap.find(new_branch_pair);
                if (it != p_branchMap.end()) {
                  printf("Found multiple lines with switched buses 1: %d 2: %d\n",
                          branch_fromBus,branch_toBus);
                  branch_index = it->second;
                  p_branchData[branch_index]->getValue(BRANCH_NUM_ELEMENTS,&nElems);
                  switched = true;
                } else {
                  boost::shared_ptr<gridpack::component::DataCollection>
                    data(new gridpack::component::DataCollection);
                  branch_index = p_branchData.size();
                  p_branchData.push_back(data);
                  nElems = 0;
                  p_branchData[branch_index]->addValue(BRANCH_NUM_ELEMENTS,nElems);
                }
            }

            if (branchAttr.first == "TransmissionElements") {
                BOOST_FOREACH( boost::property_tree::ptree::value_type lineAttr, branchAttr.second.get_child("Line"))
                {
                    loadCollection(data, lineAttr, nElems);
                }
            }
        }
    } /* end of branchAttr BOOST_FOREACH */


}

void readBus(boost::property_tree::ptree::value_type const& busesAttr, int index)
{
    boost::shared_ptr<gridpack::component::DataCollection>  data;
    int                      o_idx      = 0;
    int                      nLoads     = 0;
    int                      nGenerators = 0;
    BOOST_FOREACH( boost::property_tree::ptree::value_type busAttr, busesAttr.second)
    {
        if (busAttr.first == "Generators") {
            BOOST_FOREACH( boost::property_tree::ptree::value_type genAttr, busAttr.second.get_child("Generator"))
            {
                loadCollection(data, genAttr, nGenerators);
                ++nGenerators;
                // TODO: the current version of GOSS does not include LOAD_ID
                data->addValue(LOAD_ID, nLoads);
            }
        } else if (busAttr.first == "Loads") {
            BOOST_FOREACH( boost::property_tree::ptree::value_type loadAttr, busAttr.second.get_child("Load"))
            {
                loadCollection(data, loadAttr,nLoads);
                ++nLoads
            }
        } else {
            loadCollection(data, busAttr);
        }
        ++index;

    }

    // there is no mrid for buses. Why?

    // this changed from the original code, the insertion is done after all the
    // bus data is read and the map element is generated and inserted by retrieving
    // the BUS_NUMBER just read
    data->getValue(BUS_NUMBER,&o_idx);
    p_busData.push_back(data);
    p_busMap.insert(std::pair<int,int>(o_idx,index));

}
// ----------------------------------------------------------------
// END XML CODE
// ----------------------------------------------------------------



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
#if 0
        // debug
        printf("p[%d] BUS data size: %d\n",me,p_busData.size());
        for (i=0; i<p_busData.size(); i++) {
          printf("p[%d] Dumping bus: %d\n",me,i);
          p_busData[i]->dump();
        }
        printf("p[%d] BRANCH data size: %d\n",me,p_branchData.size());
        for (i=0; i<p_branchData.size(); i++) {
          printf("p[%d] Dumping branch: %d\n",me,i);
          p_branchData[i]->dump();
        }
#endif
        p_timer->stop(t_brdcst);
      }

private:
      /**
       * Test to see if string terminates a section
       * @return: false if first non-blank character is TERM_CHAR
       */
      bool test_end(std::string &str) const
      {
#if 1
        if (str[0] == TERM_CHAR) {
          return false;
        }
        int len = str.length();
        int i=0;
        while (i<len && str[i] == ' ') {
          i++;
        }
        if (i<len && str[i] != TERM_CHAR) {
          return true;
        } else if (i == len) {
          return true;
        } else if (str[i] == TERM_CHAR) {
          i++;
          if (i>=len || str[i] == ' ' || str[i] == '\\') {
            return false;
          } else {
            return true;
          }
        } else {
          return true;
        }
#else
        if (str[0] == '0') {
          return false;
        } else {
          return true;
        }
#endif
      }

      /**
       * Test to see if string is a comment line. Check to see if first
       * non-blank characters are "//"
       */
      bool check_comment(std::string &str) const
      {
        int ntok = str.find_first_not_of(' ',0);
        if (ntok != std::string::npos && ntok+1 != std::string::npos &&
            str[ntok] == '/' && str[ntok+1] == '/') {
          return true;
        } else {
          return false;
        }
      }
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

      // Vector of bus data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
      // Vector of branch data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
      // Map of PTI indices to index in p_busData
#ifdef OLD_MAP
      std::map<int,int> p_busMap;
#else
      boost::unordered_map<int, int> p_busMap;
#endif
      // Map of PTI index pair to index in p_branchData
#ifdef OLD_MAP
      std::map<std::pair<int, int>, int> p_branchMap;
#else
      boost::unordered_map<std::pair<int, int>, int> p_branchMap;
#endif


      std::map<std::string, XML_TYPE> typeMap; // XML
      // ----------------------------------------------------------------
      // END XML CODE
      // ----------------------------------------------------------------


      // Global variables that apply to whole network
      int p_case_id;
      double p_case_sbase;
      gridpack::utility::CoarseTimer *p_timer;
  };

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */


