/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   GOSS_test.cpp
 * @author Kevin Glass
 * @date   2014-02-11 10:28:46 d3g096
 * 
 * @brief  Test GOSS_parser capability. Currently not implemented.
 * 
 * 
 */

#include <iostream>
#include <cstdio>
#include <string>

//#include "gridpack/configuration/configuration.hpp"
//#include "gridpack/timer/coarse_timer.hpp"

#include "gridpack/parser/GOSSParser.hpp"


#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>
#include <boost/property_tree/xml_parser.hpp>

#define EPSILON     0.00001


void dumpLoads(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{

    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;

    (*data)->getValue(LOAD_BUSNUMBER, &integerValue, subIndex);
    std::cout << "LOAD_BUSNUMBER  recovered is " << integerValue << std::endl;

    (*data)->getValue(LOAD_STATUS, &integerValue, subIndex);
    std::cout << "LOAD_STATUS is " << std::endl;

    (*data)->getValue(LOAD_AREA, &integerValue, subIndex);
    std::cout << "LOAD_AREA is " << integerValue << std::endl;

    (*data)->getValue(LOAD_ZONE, &integerValue, subIndex);
    std::cout << "LOAD_ZONE is " << integerValue << std::endl;

    (*data)->getValue(LOAD_OWNER, &integerValue, subIndex);
    std::cout << "LOAD_OWNER is " << integerValue << std::endl;

    (*data)->getValue(LOAD_PL, &doubleValue, subIndex);
    std::cout << "LOAD_PL is " << doubleValue << std::endl;

    (*data)->getValue(LOAD_QL, &doubleValue, subIndex);
        std::cout << "LOAD_QL is " << doubleValue << std::endl;

    (*data)->getValue(LOAD_IP, &doubleValue, subIndex);
        std::cout << "LOAD_IP is " << doubleValue << std::endl;

    (*data)->getValue(LOAD_IQ, &doubleValue, subIndex);
        std::cout << "LOAD_IQ is " << doubleValue << std::endl;

        std::cout << "LOAD_PL is " << doubleValue << std::endl;

    (*data)->getValue(LOAD_YQ, &doubleValue, subIndex);
    std::cout << "LOAD_PL is " << doubleValue << std::endl;

}

void compareGenerators(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{

    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              stringRead;
    std::string              sourceValue;
    double                   step           = 0.0;

    (*data)->getValue(GENERATOR_BUSNUMBER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1 + subIndex*100,
           "GENERATOR_BUSNUMBER actual = " << index + 1 << " recovered is " << integerValue);

    bool test = (*data)->getValue(GENERATOR_ID, &sourceValue, subIndex);
    sprintf(stringValue, "%d", index + 1);
    stringRead  = std::string(stringValue);
    BOOST_CHECK_MESSAGE(stringRead == sourceValue,
           "GENERATOR_ID actual =" << stringRead << "-recovered is " << sourceValue);

    (*data)->getValue(GENERATOR_OWNER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_OWNER actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(GENERATOR_STAT, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_STAT actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(GENERATOR_IREG, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_IREG actual = " << index + 1 << " recovered is " << integerValue);

    step            = 1.02;
    (*data)->getValue(GENERATOR_PG, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.03;
    (*data)->getValue(GENERATOR_QG, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.04;
    (*data)->getValue(GENERATOR_QMAX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QMAX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.05;
    (*data)->getValue(GENERATOR_QMIN, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QMIN actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.06;
    (*data)->getValue(GENERATOR_VS, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_VS actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.07;
    (*data)->getValue(GENERATOR_MBASE, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_MBASE actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.08;
    (*data)->getValue(GENERATOR_ZR, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_ZR actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.09;
    (*data)->getValue(GENERATOR_ZX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_ZX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.10;
    (*data)->getValue(GENERATOR_RT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_RT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.11;
    (*data)->getValue(GENERATOR_XT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_XT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.12;
    (*data)->getValue(GENERATOR_GTAP, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_GTAP actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.13;
    (*data)->getValue(GENERATOR_RMPCT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_RMPCT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.14;
    (*data)->getValue(GENERATOR_PMAX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PMAX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.15;
    (*data)->getValue(GENERATOR_PMIN, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PMIN actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

}


BOOST_AUTO_TEST_SUITE(Parser)

void dumpBus(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{
    double                   step           = 0.0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;

    (*data)->getValue(BUS_AREA, &integerValue);
    std::cout << "BUS_AREA is " << integerValue << std::endl;

    (*data)->getValue(BUS_NUMBER, &integerValue);
    std::cout << "BUS_NUMBER is " << integerValue << std::endl;

    (*data)->getValue(BUS_TYPE, &integerValue);
    std::cout << "BUS_TYPE is " << integerValue << std::endl;

    (*data)->getValue(BUS_OWNER, &integerValue);
    std::cout << "BUS_OWNER is " << integerValue << std::endl;

    (*data)->getValue(BUS_ZONE, &integerValue);
    std::cout << "BUS_ZONE is " << integerValue << std::endl;

    (*data)->getValue(BUS_BASEKV, &doubleValue);
    std::cout << "BRANCH_FLOW_P is " << doubleValue << std::endl;

    (*data)->getValue(BUS_LOAD_PL, &doubleValue);
    std::cout << "BUS_LOAD_PL is " << doubleValue << std::endl;

    (*data)->getValue(BUS_LOAD_QL, &doubleValue);
    std::cout << "BUS_LOAD_QL is " << doubleValue << std::endl;

    (*data)->getValue(BUS_SHUNT_BL, &doubleValue);
    std::cout << "BUS_SHUNT_BL is " << doubleValue << std::endl;

    (*data)->getValue(BUS_SHUNT_GL, &doubleValue);
    std::cout << "BUS_SHUNT_GL is " << doubleValue << std::endl;

    (*data)->getValue(BUS_VOLTAGE_ANG, &doubleValue);
    std::cout << "BUS_VOLTAGE_ANG is " << doubleValue << std::endl;

    (*data)->getValue(BUS_VOLTAGE_MAG, &doubleValue);
    std::cout << "BUS_VOLTAGE_MAG is " << doubleValue << std::endl;
}


void dumpBuses(gridpack::parser::GOSSParser & parser,
    std::set<int> & hasGenerator, std::set<int> & hasLoad)
{
    std::vector <boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
#ifdef OLD_MAP
      std::map<int,int> p_busMap;
#else
      boost::unordered_map<int, int> p_busMap;
#endif
    parser.loadBuses(p_busData, p_busMap);
    int                      index          = 0;
    int                      subIndex       = 0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;
    std::set<int>::iterator  it;
    for (std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator data = p_busData.begin();
            data != p_busData.end(); ++data)
    {
        subIndex = 0;

        dumpBus(data, index, subIndex);

        while ((*data)->getValue(LOAD_BUSNUMBER, &integerValue, subIndex))
        {
//            compareLoads(data, index, subIndex);
            ++subIndex;
        }

        subIndex        = 0;
        while ((*data)->getValue(GENERATOR_BUSNUMBER, &integerValue, subIndex))
        {
//            compareGenerators(data, index, subIndex);
            ++subIndex;
        }
        ++index;
    }
}


void dumpBranches(gridpack::parser::GOSSParser & parser)
{
    int                      index          = 0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              stringRead;
    std::string              sourceValue;
    double                   step           = 0.0;

    // load branches from parser
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > branches;
    parser.loadBranches(branches);

    for (std::vector<boost::shared_ptr<gridpack::component::DataCollection> >::iterator data = branches.begin();
            data < branches.end(); data++)
    {

        (*data)->getValue(BRANCH_FROMBUS, &integerValue);
        std::cout << "BRANCH_FROMBUS is " << integerValue << std::endl;

        (*data)->getValue(BRANCH_TOBUS, &integerValue);
        std::cout << "BRANCH_TOBUS is " << integerValue << std::endl;

        (*data)->getValue(BRANCH_INDEX, &integerValue,0);
        std::cout << "BRANCH_INDEX is " << integerValue << std::endl;

        (*data)->getValue("mrid", &sourceValue, 0);
        sprintf(stringValue, "%d", index + 1);
        stringRead  = std::string(stringValue);
        std::cout << "mrid is " << sourceValue << std::endl;

        step            = 1.01;
        (*data)->getValue(BRANCH_FLOW_P, &doubleValue,0);
            std::cout << "BRANCH_FLOW_P is " << doubleValue << std::endl;

        step            = 1.02;
        (*data)->getValue(BRANCH_FLOW_Q, &doubleValue,0);
            std::cout << "BRANCH_FLOW_Q is " << doubleValue << std::endl;

        step            = 1.03;
        (*data)->getValue(BRANCH_R, &doubleValue,0);
            std::cout << "BRANCH_R is " << doubleValue << std::endl;

        step            = 1.04;
        (*data)->getValue(BRANCH_RATING_A, &doubleValue,0);
            std::cout << "BRANCH_RATING_A is " << doubleValue << std::endl;

        step            = 1.05;
        (*data)->getValue(BRANCH_RATING_B, &doubleValue,0);
            std::cout << "BRANCH_RATING_B is " << doubleValue << std::endl;

        step            = 1.06;
        (*data)->getValue(BRANCH_RATING_C, &doubleValue,0);
            std::cout << "BRANCH_RATING_C is " << doubleValue << std::endl;

        step            = 1.07;
        (*data)->getValue(BRANCH_RATING, &doubleValue,0);
            std::cout << "BRANCH_RATING is " << doubleValue << std::endl;

        (*data)->getValue(BRANCH_STATUS, &integerValue,0);
        std::cout << "BRANCH_STATUS is " << integerValue << std::endl;

        step            = 1.08;
        (*data)->getValue(BRANCH_X, &doubleValue,0);
            std::cout << "BRANCH_X is " << doubleValue << std::endl;

        step            = 1.09;
        (*data)->getValue(BRANCH_B, &doubleValue,0);
            std::cout << "BRANCH_B is " << doubleValue << std::endl;

        step            = 1.10;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_B1, &doubleValue,0);
            std::cout << "BRANCH_SHUNT_ADMTTNC_B1 is " << doubleValue << std::endl;

        step            = 1.11;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_B2, &doubleValue,0);
            std::cout << "BRANCH_SHUNT_ADMTTNC_B2 is " << doubleValue << std::endl;

        step            = 1.12;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_G1, &doubleValue,0);
            std::cout << "BRANCH_SHUNT_ADMTTNC_G1 is " << doubleValue << std::endl;

        step            = 1.13;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_G2, &doubleValue,0);
            std::cout << "BRANCH_SHUNT_ADMTTNC_G2 is " << doubleValue << std::endl;

        ++index;
    }
}
void compareLoads(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{

    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;

    (*data)->getValue(LOAD_BUSNUMBER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1 + 100*subIndex,
           "LOAD_BUSNUMBER actual = " << index + 1 + 100*subIndex << " recovered is " << integerValue);

    (*data)->getValue(LOAD_STATUS, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "LOAD_STATUS actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(LOAD_AREA, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "LOAD_AREA actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(LOAD_ZONE, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "LOAD_ZONE actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(LOAD_OWNER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "LOAD_OWNER actual = " << index + 1 << " recovered is " << integerValue);

    step            = 1.18;
    (*data)->getValue(LOAD_PL, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_PL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.19;
    (*data)->getValue(LOAD_QL, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_QL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.20;
    (*data)->getValue(LOAD_IP, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_IP actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.21;
    (*data)->getValue(LOAD_IQ, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_IQ actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.22;
    (*data)->getValue(LOAD_YP, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_PL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.23;
    (*data)->getValue(LOAD_YQ, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "LOAD_PL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

}

void compareGenerators(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{

    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              stringRead;
    std::string              sourceValue;
    double                   step           = 0.0;

    (*data)->getValue(GENERATOR_BUSNUMBER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1 + subIndex*100,
           "GENERATOR_BUSNUMBER actual = " << index + 1 << " recovered is " << integerValue);

    bool test = (*data)->getValue(GENERATOR_ID, &sourceValue, subIndex);
    sprintf(stringValue, "%d", index + 1);
    stringRead  = std::string(stringValue);
    BOOST_CHECK_MESSAGE(stringRead == sourceValue,
           "GENERATOR_ID actual =" << stringRead << "-recovered is " << sourceValue);

    (*data)->getValue(GENERATOR_OWNER, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_OWNER actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(GENERATOR_STAT, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_STAT actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(GENERATOR_IREG, &integerValue, subIndex);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "GENERATOR_IREG actual = " << index + 1 << " recovered is " << integerValue);

    step            = 1.02;
    (*data)->getValue(GENERATOR_PG, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.03;
    (*data)->getValue(GENERATOR_QG, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.04;
    (*data)->getValue(GENERATOR_QMAX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QMAX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.05;
    (*data)->getValue(GENERATOR_QMIN, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_QMIN actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.06;
    (*data)->getValue(GENERATOR_VS, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_VS actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.07;
    (*data)->getValue(GENERATOR_MBASE, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_MBASE actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.08;
    (*data)->getValue(GENERATOR_ZR, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_ZR actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.09;
    (*data)->getValue(GENERATOR_ZX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_ZX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.10;
    (*data)->getValue(GENERATOR_RT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_RT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.11;
    (*data)->getValue(GENERATOR_XT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_XT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.12;
    (*data)->getValue(GENERATOR_GTAP, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_GTAP actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.13;
    (*data)->getValue(GENERATOR_RMPCT, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_RMPCT actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.14;
    (*data)->getValue(GENERATOR_PMAX, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PMAX actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.15;
    (*data)->getValue(GENERATOR_PMIN, &doubleValue, subIndex);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "GENERATOR_PMIN actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

}

void compareBus(std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator & data,
    int index, int subIndex)
{
    double                   step           = 0.0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;

    (*data)->getValue(BUS_AREA, &integerValue);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "BUS_AREA actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(BUS_NUMBER, &integerValue);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "BUS_NUMBER actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(BUS_TYPE, &integerValue);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "BUS_TYPE actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(BUS_OWNER, &integerValue);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "BUS_OWNER actual = " << index + 1 << " recovered is " << integerValue);

    (*data)->getValue(BUS_ZONE, &integerValue);
    BOOST_CHECK_MESSAGE(integerValue == index + 1,
           "BUS_ZONE actual = " << index + 1 << " recovered is " << integerValue);

    step            = 1.01;
    (*data)->getValue(BUS_BASEKV, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BRANCH_FLOW_P actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.16;
    (*data)->getValue(BUS_LOAD_PL, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_LOAD_PL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.17;
    (*data)->getValue(BUS_LOAD_QL, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_LOAD_QL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.24;
    (*data)->getValue(BUS_SHUNT_BL, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_SHUNT_BL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.25;
    (*data)->getValue(BUS_SHUNT_GL, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_SHUNT_GL actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.26;
    (*data)->getValue(BUS_VOLTAGE_ANG, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_VOLTAGE_ANG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;

    step            = 1.27;
    (*data)->getValue(BUS_VOLTAGE_MAG, &doubleValue);
    BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
    if ((double)index + step > doubleValue + EPSILON ||
            (double)index + step < doubleValue - EPSILON )
        std::cout << "BUS_VOLTAGE_MAG actual = " << (double)index + step << " recovered is " << doubleValue << std::endl;
}

void compareBuses(gridpack::parser::GOSSParser & parser,
    std::set<int> & hasGenerator, std::set<int> & hasLoad)
{
    std::vector <boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
#ifdef OLD_MAP
      std::map<int,int> p_busMap;
#else
      boost::unordered_map<int, int> p_busMap;
#endif
    parser.loadBuses(p_busData, p_busMap);
    int                      index          = 0;
    int                      subIndex       = 0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;
    std::set<int>::iterator  it;
    for (std::vector <boost::shared_ptr<gridpack::component::DataCollection > >::iterator data = p_busData.begin();
            data != p_busData.end(); ++data)
    {
        subIndex = 0;

        dumpBus(data, index, subIndex);

        while ((*data)->getValue(LOAD_BUSNUMBER, &integerValue, subIndex))
        {
//            compareLoads(data, index, subIndex);
            ++subIndex;
        }

        subIndex        = 0;
        while ((*data)->getValue(GENERATOR_BUSNUMBER, &integerValue, subIndex))
        {
//            compareGenerators(data, index, subIndex);
            ++subIndex;
        }
        ++index;
    }
}

void compareBranches(gridpack::parser::GOSSParser & parser)
{
    int                      index          = 0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              stringRead;
    std::string              sourceValue;
    double                   step           = 0.0;

    // load branches from parser
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > branches;
    parser.loadBranches(branches);

    for (std::vector<boost::shared_ptr<gridpack::component::DataCollection> >::iterator data = branches.begin();
            data < branches.end(); data++)
    {

        (*data)->getValue(BRANCH_FROMBUS, &integerValue);
        BOOST_CHECK_MESSAGE(integerValue == 2*index + 1,
               "BRANCH_FROMBUS actual = " << 2*index + 1 << " recovered is " << integerValue);

        (*data)->getValue(BRANCH_TOBUS, &integerValue);
        BOOST_CHECK_MESSAGE(integerValue == 2*index + 2,
               "BRANCH_TOBUS actual = " << 2*index + 2 << " recovered is " << integerValue);

        (*data)->getValue(BRANCH_INDEX, &integerValue,0);
        BOOST_CHECK_MESSAGE(integerValue == index + 1,
               "BRANCH_INDEX actual = " << index + 1 << " recovered is " << integerValue);

        (*data)->getValue("mrid", &sourceValue, 0);
        sprintf(stringValue, "%d", index + 1);
        stringRead  = std::string(stringValue);
        BOOST_CHECK_MESSAGE(stringRead == sourceValue,
            "mrid actual = " << stringRead << " recovered is " << sourceValue);

        step            = 1.01;
        (*data)->getValue(BRANCH_FLOW_P, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_FLOW_P actual = " << (double)index + 1.01 << " recovered is " << doubleValue << std::endl;

        step            = 1.02;
        (*data)->getValue(BRANCH_FLOW_Q, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_FLOW_Q actual = " << (double)index + 1.02 << " recovered is " << doubleValue << std::endl;

        step            = 1.03;
        (*data)->getValue(BRANCH_R, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_R actual = " << (double)index + 1.03 << " recovered is " << doubleValue << std::endl;

        step            = 1.04;
        (*data)->getValue(BRANCH_RATING_A, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_RATING_A actual = " << (double)index + 1.04 << " recovered is " << doubleValue << std::endl;

        step            = 1.05;
        (*data)->getValue(BRANCH_RATING_B, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_RATING_B actual = " << (double)index + 1.05 << " recovered is " << doubleValue << std::endl;

        step            = 1.06;
        (*data)->getValue(BRANCH_RATING_C, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_RATING_C actual = " << (double)index + 1.06 << " recovered is " << doubleValue << std::endl;

        step            = 1.07;
        (*data)->getValue(BRANCH_RATING, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_RATING actual = " << (double)index + 1.07 << " recovered is " << doubleValue << std::endl;

        (*data)->getValue(BRANCH_STATUS, &integerValue,0);
        BOOST_CHECK_MESSAGE(integerValue == index + 1,
            "BRANCH_STATUS actual = " << index +1 << " recovered is " << integerValue);

        step            = 1.08;
        (*data)->getValue(BRANCH_X, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_X actual = " << (double)index + 1.08 << " recovered is " << doubleValue << std::endl;

        step            = 1.09;
        (*data)->getValue(BRANCH_B, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_B actual = " << (double)index + 1.09 << " recovered is " << doubleValue << std::endl;

        step            = 1.10;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_B1, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_SHUNT_ADMTTNC_B1 actual = " << (double)index + 1.1 << " recovered is " << doubleValue << std::endl;

        step            = 1.11;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_B2, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_SHUNT_ADMTTNC_B2 actual = " << (double)index + 1.11 << " recovered is " << doubleValue << std::endl;

        step            = 1.12;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_G1, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_SHUNT_ADMTTNC_G1 actual = " << (double)index + 1.12 << " recovered is " << doubleValue << std::endl;

        step            = 1.13;
        (*data)->getValue(BRANCH_SHUNT_ADMTTNC_G2, &doubleValue,0);
        BOOST_CHECK_CLOSE(doubleValue, (double)index + step, EPSILON);
        if ((double)index + step > doubleValue + EPSILON ||
                (double)index + step < doubleValue - EPSILON )
            std::cout << "BRANCH_SHUNT_ADMTTNC_G2 actual = " << (double)index + 1.13 << " recovered is " << doubleValue << std::endl;

        ++index;
    }

}

BOOST_AUTO_TEST_CASE(ArtificialData)
{

    bool ok = true;
    std::string   fileName("gridpack-test1.xml");
    std::cout << "TESTING LOAD OF " << fileName << std::endl;
    gridpack::parser::GOSSParser           parser;

    try {
        parser.parse(fileName.c_str());
    } catch (boost::property_tree::xml_parser_error & e) {
        e.what();
    } catch (boost::exception & e) {
        std::cout << "General exception\n\t";
        std::cout << std::endl;
    }

    int busesWithGenerators[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int busesWithLoad[]       = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::set<int> hasGenerator(busesWithGenerators, busesWithGenerators+10);
    std::set<int> hasLoad(busesWithLoad, busesWithLoad+10);


//    dumpBuses(parser, hasGenerator, hasLoad);
    compareBuses(parser, hasGenerator, hasLoad);
//    dumpBranches(parser);
    compareBranches(parser);
    std::cout << "END TESTING LOAD OF " << fileName << std::endl;


}

BOOST_AUTO_TEST_CASE(ArtificialData2)
{
    bool ok = true;
    std::string   fileName("gridpack-test2.xml");
    std::cout << "TESTING LOAD OF " << fileName << std::endl;
    gridpack::parser::GOSSParser           parser;

    try {
        parser.parse(fileName.c_str());
    } catch (boost::property_tree::xml_parser_error & e) {
        e.what();
    } catch (boost::exception & e) {
        std::cout << "General excpetion\n\t";
    }
//    std::string caseId      = parser.getCaseId();
//    std::string caseSbase("");

    BOOST_CHECK_MESSAGE(parser.getCaseId() == std::string("Case 1"),
           "BRANCH_FROMBUS actual = Case 1 recovered is " << parser.getCaseId());

    int busesWithGenerators[] = {0, 1, 2, 3, 4, 5, 6, 8, 9};
    int busesWithLoad[]       = {0, 1, 2, 3, 4, 5, 7, 8, 9};
    std::set<int> hasGenerator(busesWithGenerators, busesWithGenerators+9);
    std::set<int> hasLoad(busesWithLoad, busesWithLoad+10);
    compareBuses(parser, hasGenerator, hasLoad);
    compareBranches(parser);
    std::cout << "END TESTING LOAD OF " << fileName << std::endl;
}

BOOST_AUTO_TEST_CASE(Greek_118)
{

    std::string   fileName("full_gp_north_118.xml");
    std::cout << "TESTING LOAD OF " << fileName << std::endl;
    gridpack::parser::GOSSParser           parser;

    try {
        parser.parse(fileName.c_str());
    } catch (boost::property_tree::xml_parser_error & e) {
        e.what();
    } catch (boost::exception & e) {
        std::cout << "General excpetion\n\t";
        std::cout << std::endl;
    }
//    dumpBranches(parser);
    std::cout << "END TESTING LOAD OF " << fileName << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{

    return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
    int                      result         = 0;
    printf("Testing Aritificial Input\n");
    result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );

    return result;
}


