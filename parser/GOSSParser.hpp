/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * GOSSparser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef GOSSPARSER_HPP_
#define GOSSPARSER_HPP_

#define OLD_MAP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include <vector>

#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace parser {

typedef std::vector <std::vector <boost::shared_ptr<gridpack::component::DataCollection > > > BucketList;

class GOSSParser
{
    public:

    enum XML_TYPE {INTEGER, BOOLEAN, DOUBLE, CHARACTER, STRING};

    /// Constructor
    explicit GOSSParser() : busIterator(0), branchIterator(0),  nBuses(0), nBranches(0), maxBus(0), p_case_id(1), p_case_sbase(1.0){};


    /**
     * Destructor
     */
    virtual ~GOSSParser(){}

    /**
     * Parse network configuration file and create network
     * @param fileName name of network file
     */
    virtual void parse(const std::string & xmlFileName)
    {
        // setup parser
        boost::property_tree::ptree xmlTree;

        try {
            setupXMLTree(xmlFileName, xmlTree);
        } catch (boost::property_tree::xml_parser_error & e) {
            std::cout << "ERROR IN PARSER SETUP" << std::endl;
            throw e;
        }

        setTypeAssociations(xmlTree);
        setCase();
        loadBusData(xmlTree);
        loadBranchData(xmlTree);
    }

    void loadBranches(std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & branches)
    {
        for (int i = 0; i < maxBus; ++i)
        {
            for (int j = 0; j < maxBus; ++j)
            {
                if (branchBucket[i][j]) {
                    branches.push_back(branchBucket[i][j]);
                }
            }
        }
    }
#ifdef OLD_MAP
    void loadBuses(std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & buses,
        std::map<int,int> busMap)
#else
    void loadBuses(std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & buses,
      boost::unordered_map<int, int> busMap)
#endif
    {
        buses = p_busData;
        busMap = p_busMap;
    }

    int getNBuses(){return nBuses;};
    int getNBranches(){return nBranches;};
    int getCaseId(){return p_case_id;};
    int getCaseSbase(){return p_case_sbase;};

    protected:
    private:
    // ****************************************************************************
    // ****************************************************************************

    void setupXMLTree(const std::string & xmlFileName,
        boost::property_tree::ptree & xmlTree)
    {
        try {
            read_xml(xmlFileName, xmlTree);
        } catch (boost::property_tree::xml_parser_error  & e) {
            throw e;
        }
    }

    void setCase()
    {
        p_case_id       = 1;
        p_case_sbase    = 1.0;
    }


    void setTypeAssociations(boost::property_tree::ptree &  xmlTree)
    {
        std::string          name("");
        boost::property_tree::ptree schemaStyleAttr;

        try {
            schemaStyleAttr =
                    xmlTree.get_child("application.grammars.xs:schema");
        } catch (boost::exception & e){
            std::cout << "Throwing exception from get_child(\"application.grammars.xs:schema\")" << std::endl;
            throw;
        }

        typeMap["mrid"] = STRING;
        typeMap["elementIndex"] = INTEGER;

        BOOST_FOREACH( boost::property_tree::ptree::value_type const& elementAttr, schemaStyleAttr)
        {
            try {
                recursiveSetTypeAssociation(elementAttr, name);
            } catch (boost::exception & e) {
                std::cout << "Error in recursive type setting" << std::endl;
            }
        }
    }

    void recursiveSetTypeAssociation(boost::property_tree::ptree::value_type const& elementAttr, std::string & name)
    {
        XML_TYPE             typeEnum;
        // if the XML type starts with xs:element, then it is component value
        if (elementAttr.first == "xs:element") {
            const boost::property_tree::ptree & attributes = elementAttr.second.get_child("<xmlattr>");
            BOOST_FOREACH( boost::property_tree::ptree::value_type const& attr, attributes)
            {
                std::string               type   = attr.second.data();
                std::vector<std::string>  split_type;
                boost::algorithm::split(split_type, type, boost::algorithm::is_any_of(":"), boost::token_compress_on);
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
                recursiveSetTypeAssociation(sequenceAttr, name);
            }
        }
    }

    void loadBusData(boost::property_tree::ptree & xmlTree)
    {
        boost::property_tree::ptree busXMLTree =
                xmlTree.get_child("application.GridpackPowergrid.Buses");
        int                  nBus           = 0;
        BOOST_FOREACH( boost::property_tree::ptree::value_type const& busesAttr, busXMLTree)
        {
            readBus(busesAttr);
            ++nBus;
        }
    }

    void readBus(boost::property_tree::ptree::value_type const & busesAttr)
    {
        boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);

        int                      nLoads     = 0;
        int                      nGenerators = 0;

        BOOST_FOREACH( boost::property_tree::ptree::value_type busAttr,
            busesAttr.second)
        {
            if (busAttr.first == "Generators")
            {
                BOOST_FOREACH( boost::property_tree::ptree::value_type genSet,
                        busAttr.second)
                {
                    BOOST_FOREACH( boost::property_tree::ptree::value_type genAttr,
                            genSet.second)
                    {
                          loadCollection(data, genAttr, nGenerators);
                    }
                    ++nGenerators;
                }

            } else if (busAttr.first == "Loads") {

                BOOST_FOREACH( boost::property_tree::ptree::value_type loadSet,
                        busAttr.second)
                {
                    BOOST_FOREACH( boost::property_tree::ptree::value_type loadAttr,
                            loadSet.second)
                    {
                        loadCollection(data, loadAttr, nLoads);
                    }
                    ++nLoads;
                }

            } else {
                loadCollection(data, busAttr);
                int      busIndex = 0;

                data->getValue(BUS_NUMBER, &busIndex);
                if (maxBus < busIndex) maxBus = busIndex;
                p_busMap.insert(std::pair<int,int>(busIndex,busIndex));
            }
        }
        p_busData.push_back(data);
        ++nBuses;
    }

    void loadBranchData(boost::property_tree::ptree & xmlTree)
    {
        boost::property_tree::ptree branchesXMLTree;
        branchesXMLTree =
                xmlTree.get_child("application.GridpackPowergrid.Branches");
        int                  index          = 0;

        // create an nBuses x nBuses bucket list
        ++maxBus;
        branchBucket    = BucketList(maxBus, std::vector<boost::shared_ptr<gridpack::component::DataCollection> >(maxBus));

        BOOST_FOREACH( boost::property_tree::ptree::value_type const& branchesAttr, branchesXMLTree)
        {
            readBranch(branchesAttr, index);
            ++ index;
        }
    }

    void readBranch(boost::property_tree::ptree::value_type const& branchesAttr, int index)
    {
        int                      nLoads     = 0;
        int                      nGenerators = 0;
        std::map<std::string, std::string>     keyValue;
        int              nElems       = 0;

        // load the branch data into the keyValue map
        BOOST_FOREACH( boost::property_tree::ptree::value_type branchAttr, branchesAttr.second)
        {
            try{

                if (branchAttr.first == "TransmissionElements") {
                    // read xml branch and store in map
                    BOOST_FOREACH( boost::property_tree::ptree::value_type lineAttr, branchAttr.second.get_child("Line"))
                    {
                            keyValue[lineAttr.first] = lineAttr.second.data();
                    }
                } else {
                    if (branchAttr.first == BRANCH_FROMBUS || branchAttr.first == BRANCH_TOBUS)
                    {
                        keyValue[branchAttr.first] = branchAttr.second.data();
                    };
                }
            }catch (boost::exception & e) {
                std::cout << "Failed to read transmission element" << std::endl;
                throw;
            }
        }

        // branch from/to indices
        int       fromBus = atoi(keyValue["BRANCH_FROMBUS"].c_str());
        int       toBus   = atoi(keyValue["BRANCH_TOBUS"].c_str());

        // bucketList load branch bucket using from/to or to/from for switched buses
        if (branchBucket[fromBus][toBus])
        {
            branchBucket[fromBus][toBus]->getValue(BRANCH_NUM_ELEMENTS,&nElems);
        } else if (branchBucket[toBus][fromBus]) {
            branchBucket[toBus][fromBus]->getValue(BRANCH_NUM_ELEMENTS,&nElems);
        } else {
            nElems      = 0;
            boost::shared_ptr<gridpack::component::DataCollection>  s_data(new gridpack::component::DataCollection);
            s_data->addValue(BRANCH_FROMBUS, fromBus);
            s_data->addValue(BRANCH_TOBUS, toBus);

            branchBucket[fromBus][toBus] = s_data;
        }

        // load key/value elements into bucket
        loadKeyValues(keyValue, branchBucket[fromBus][toBus], nElems);

        // set the element index in the data collection
        ++nElems;
        branchBucket[fromBus][toBus]->setValue(BRANCH_NUM_ELEMENTS,nElems);

        ++nBranches;
    }

    void loadKeyValues(std::map<std::string,std::string> &  keyValue,
        boost::shared_ptr<gridpack::component::DataCollection> data, int nElems)
    {
        for (std::map<std::string,std::string>::iterator it = keyValue.begin();
                it != keyValue.end(); ++it)
        {
            XML_TYPE         type               = typeMap[it->first];

            switch (type)
            {
                case BOOLEAN:
                {
                    bool             value          = true;
                    if (it->second == "false") value = false;
                    data->addValue(it->first.c_str(), value, nElems);
                    break;
                }
                case INTEGER:
                {
                    data->addValue(it->first.c_str(), atoi(it->second.c_str()), nElems);
                    break;
                }
                case DOUBLE:
                {
                    data->addValue(it->first.c_str(), atof(it->second.c_str()), nElems);
                    break;
                }
                case CHARACTER:
                {
                    data->addValue(it->first.c_str(), it->second[0], nElems);
                    break;
                }
                case STRING:
                {
                    data->addValue(it->first.c_str(), it->second.c_str(), nElems);
                    break;
                }
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

    int                      busIterator;
    int                      branchIterator;
    int                      nBuses;
    int                      nBranches;
    int                      maxBus;
    // Vector of bus data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
                             p_busData;
#ifdef OLD_MAP
      std::map<int,int> p_busMap;
#else
      boost::unordered_map<int, int> p_busMap;
#endif
    // Vector of branch data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
                             p_branchData;

    // Map of PTI indices to index in p_busData
    BucketList               branchBucket;

    int                      p_case_id;
    double                   p_case_sbase;
    std::map<std::string, XML_TYPE> typeMap;

}; /* end of GOSSParser */

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

} /* namespace parser */
} /* namespace gridpack */

#endif /* GOSSPARSER_HPP_ */
