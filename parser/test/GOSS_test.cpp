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
//#include <boost/property_tree/xml_parser.hpp>

#define EPSILON     0.00001



bool compareString(const char * str1, const char * str2)
{
    int                      len            = strlen(str1) + 2;
    bool                     comp           = true;

    std::cout << "char\t" << str1 << "\t" << str2 << "\tDIFF" << std::endl;
    for (int i = 0; i < len; i++)
    {
        std::cout << i << "\t" << (int)str1[i] << "\t" << (int)str2[i] << "\t" << (int)str1[i] - (int)str2[i] <<std::endl;
        if (str1[i] != str2[i]) comp = false;
    }

    return comp;
}

BOOST_AUTO_TEST_SUITE(Parser)

/*
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

*/
bool verifyBusData(gridpack::component::DataCollection & data,
    gridpack::parser::GOSSParser parser, int busIndex)
{
    bool                     valid          = true;
    bool                     busValid       = true;
    double                   step           = 0.0;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;

    data.getValue(BUS_AREA, &integerValue);
    valid               = integerValue == busIndex + 1;
    busValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BUS_AREA actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(BUS_NUMBER, &integerValue);
    valid               = integerValue == busIndex + 1;
    busValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BUS_NUMBER actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(BUS_TYPE, &integerValue);
    valid               = integerValue == busIndex + 1;
    busValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BUS_TYPE actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(BUS_OWNER, &integerValue);
    valid               = integerValue == busIndex + 1;
    busValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BUS_OWNER actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(BUS_ZONE, &integerValue);
    valid               = integerValue == busIndex + 1;
    busValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BUS_ZONE actual = " << busIndex + 1 << " recovered is " << integerValue);

    step            = 1.01;
    data.getValue(BUS_BASEKV, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_BASEKV actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.16;
    data.getValue(BUS_LOAD_PL, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_LOAD_PL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.17;
    data.getValue(BUS_LOAD_QL, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_LOAD_QL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.24;
    data.getValue(BUS_SHUNT_BL, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_SHUNT_BL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.25;
    data.getValue(BUS_SHUNT_GL, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_SHUNT_GL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.26;
    data.getValue(BUS_VOLTAGE_ANG, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_VOLTAGE_ANG actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.27;
    data.getValue(BUS_VOLTAGE_MAG, &doubleValue);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid           |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BUS_VOLTAGE_MAG actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    return busValid;
}

bool verifyGeneratorData(gridpack::component::DataCollection & data,
    gridpack::parser::GOSSParser parser, int busIndex, int generatorIndex)
{
    bool                     valid          = true;
    bool                     busValid       = true;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    char                     stringValue[64];
    std::string              stringRead;
    std::string              sourceValue;
    double                   step           = 0.0;

    data.getValue(GENERATOR_BUSNUMBER, &integerValue, generatorIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "GENERATOR_BUSNUMBER actual = " << busIndex + 1 << " recovered is " << integerValue);

    bool test = data.getValue(GENERATOR_ID, &sourceValue, generatorIndex);
    sprintf(stringValue, "%d", busIndex + 1);
    valid               = !strcmp(stringValue, sourceValue.c_str());
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "GENERATOR_ID actual = " << stringValue << " recovered is " << sourceValue);

    data.getValue(GENERATOR_OWNER, &integerValue, generatorIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "GENERATOR_OWNER actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(GENERATOR_STAT, &integerValue, generatorIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "GENERATOR_STAT actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(GENERATOR_IREG, &integerValue, generatorIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "GENERATOR_IREG actual = " << busIndex + 1 << " recovered is " << integerValue);

    step            = 1.02;
    data.getValue(GENERATOR_PG, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_PG actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.03;
    data.getValue(GENERATOR_QG, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_QG actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.04;
    data.getValue(GENERATOR_QMAX, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_QMAX actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.05;
    data.getValue(GENERATOR_QMIN, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_QMIN actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.06;
    data.getValue(GENERATOR_VS, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_VS actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.07;
    data.getValue(GENERATOR_MBASE, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_MBASE actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.08;
    data.getValue(GENERATOR_ZR, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_ZR actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.09;
    data.getValue(GENERATOR_ZX, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_ZX actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.10;
    data.getValue(GENERATOR_RT, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_RT actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.11;
    data.getValue(GENERATOR_XT, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_XT actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.12;
    data.getValue(GENERATOR_GTAP, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_GTAP actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.13;
    data.getValue(GENERATOR_RMPCT, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_RMPCT actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.14;
    data.getValue(GENERATOR_PMAX, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_PMAX actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.15;
    data.getValue(GENERATOR_PMIN, &doubleValue, generatorIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "GENERATOR_PMIN actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    return busValid;
};

bool verifyLoadData(gridpack::component::DataCollection & data,
    gridpack::parser::GOSSParser parser, int busIndex, int loadIndex)
{
    bool                     valid          = true;
    bool                     busValid       = true;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    bool                     logicalValue   = true;
    char                     charValue      = 'a';
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;

    data.getValue(LOAD_BUSNUMBER, &integerValue, loadIndex);
    valid               = integerValue == busIndex + 1 + 100*loadIndex;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "LOAD_BUSNUMBER actual = " << busIndex + 1 + 100*loadIndex << " recovered is " << integerValue);

    data.getValue(LOAD_STATUS, &integerValue, loadIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "LOAD_STATUS actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(LOAD_AREA, &integerValue, loadIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "LOAD_AREA actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(LOAD_ZONE, &integerValue, loadIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "LOAD_ZONE actual = " << busIndex + 1 << " recovered is " << integerValue);

    data.getValue(LOAD_OWNER, &integerValue, loadIndex);
    valid               = integerValue == busIndex + 1;
    busValid           |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "LOAD_OWNER actual = " << busIndex + 1 << " recovered is " << integerValue);

    step            = 1.18;
    data.getValue(LOAD_PL, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_PL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.19;
    data.getValue(LOAD_QL, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_QL actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.20;
    data.getValue(LOAD_IP, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_IP actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.21;
    data.getValue(LOAD_IQ, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_IQ actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.22;
    data.getValue(LOAD_YP, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_YP actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    step            = 1.23;
    data.getValue(LOAD_YQ, &doubleValue, loadIndex);
    valid               =
            ((double)(busIndex + step) < doubleValue + EPSILON &&
                    (double)(busIndex + step) > doubleValue - EPSILON );
    busValid              |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)busIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "LOAD_YQ actual = " << (double)busIndex + step << " recovered is " << doubleValue << std::endl;
    }

    return busValid;
};


bool validateBus(gridpack::component::DataCollection & data,
    gridpack::parser::GOSSParser parser, int busIndex)
{
    bool                     busValid       = true;
    int                      generatorIndex = 0;
    int                      loadIndex      = 0;
    int                      integerValue   = 0;

    busValid           |= verifyBusData(data, parser, busIndex);
    while (data.getValue(GENERATOR_BUSNUMBER, &integerValue, generatorIndex))
    {
        busValid        |= verifyGeneratorData(data, parser, busIndex, generatorIndex);
        ++generatorIndex;
    }

    while (data.getValue(LOAD_BUSNUMBER, &integerValue, loadIndex))
    {
        busValid        |= verifyLoadData(data, parser, busIndex, loadIndex);
        ++loadIndex;
    }

    return busValid;
}

bool validateBranch(gridpack::component::DataCollection & data,
    gridpack::parser::GOSSParser parser, int branchIndex)
{
    bool                     valid          = true;
    bool                     branchValid    = true;
    int                      integerValue   = 0;
    double                   doubleValue    = 0.0;
    char                     stringValue[64];
    std::string              sourceValue;
    double                   step           = 0.0;


//    parser.test_dumpDataColletion(data);
    // validate BRANCH_FROMBUS
    data.getValue(BRANCH_FROMBUS, &integerValue);
    valid               = integerValue == 2*branchIndex + 1;
    branchValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BRANCH_FROMBUS actual = " << 2*branchIndex + 1 << " recovered is " << integerValue);

    // validate BRANCH_TOBUS
    data.getValue(BRANCH_TOBUS, &integerValue);
    valid               = integerValue == 2*branchIndex + 2;
    branchValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BRANCH_TOBUS actual = " << 2*branchIndex + 2 << " recovered is " << integerValue);

    // validate BRANCH_INDEX
    data.getValue(BRANCH_INDEX, &integerValue);
    valid               = integerValue == branchIndex + 1;
    branchValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BRANCH_INDEX actual = " << branchIndex + 1 << " recovered is " << integerValue);

    // validate mrid
    data.getValue("mrid", &sourceValue);
    sprintf(stringValue, "%d", branchIndex + 1);
    valid               = !strcmp(stringValue, sourceValue.c_str());
    branchValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "mrid actual = " << stringValue << " recovered is " << sourceValue);

    // validate BRANCH_FLOW_P
    step            = 1.01;
    data.getValue(BRANCH_FLOW_P, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_FLOW_P actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_FLOW_Q
    step            = 1.02;
    data.getValue(BRANCH_FLOW_Q, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_FLOW_Q actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_R
    step            = 1.03;
    data.getValue(BRANCH_R, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_R actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_RATING_A
    step            = 1.04;
    data.getValue(BRANCH_RATING_A, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_RATING_A actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_RATING_B
    step            = 1.05;
    data.getValue(BRANCH_RATING_B, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_RATING_B actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_RATING_C
    step            = 1.06;
    data.getValue(BRANCH_RATING_C, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_RATING_C actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_RATING
    step            = 1.07;
    data.getValue(BRANCH_RATING, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_RATING actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_STATUS
    data.getValue(BRANCH_STATUS, &integerValue);
    valid               = integerValue == branchIndex + 1;
    branchValid        |= valid;
    BOOST_CHECK_MESSAGE(valid,
           "BRANCH_STATUS actual = " << branchIndex + 1 << " recovered is " << integerValue);

    // validate BRANCH_X
    step            = 1.08;
    data.getValue(BRANCH_X, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_X actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_B
    step            = 1.09;
    data.getValue(BRANCH_B, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_B actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_SHUNT_ADMTTNC_B1
    step            = 1.10;
    data.getValue(BRANCH_SHUNT_ADMTTNC_B1, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_SHUNT_ADMTTNC_B1 actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_SHUNT_ADMTTNC_B2
    step            = 1.11;
    data.getValue(BRANCH_SHUNT_ADMTTNC_B2, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_SHUNT_ADMTTNC_B2 actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_SHUNT_ADMTTNC_G1
    step            = 1.12;
    data.getValue(BRANCH_SHUNT_ADMTTNC_G1, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_SHUNT_ADMTTNC_G1 actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    // validate BRANCH_SHUNT_ADMTTNC_G2
    step            = 1.13;
    data.getValue(BRANCH_SHUNT_ADMTTNC_G2, &doubleValue);
    valid               =
            ((double)(branchIndex + step) < doubleValue + EPSILON &&
                    (double)(branchIndex + step) > doubleValue - EPSILON );
    branchValid        |= valid;
    BOOST_CHECK_CLOSE(doubleValue, (double)branchIndex + step, EPSILON);
    if(!valid)
    {
        std::cout << "BRANCH_SHUNT_ADMTTNC_G2 actual = " << (double)branchIndex + step << " recovered is " << doubleValue << std::endl;
    }

    return branchValid;
}

void validateBusCollection(std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & busCollectionVector,
    gridpack::parser::GOSSParser & parser)
{
    bool                     busValid       = true;
    int                      busIndex       = 0;

    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >::iterator dataCollection;
    for (dataCollection = busCollectionVector.begin();
            dataCollection != busCollectionVector.end(); ++dataCollection)
    {
        busValid        = validateBus(**dataCollection, parser, busIndex);
        ++busIndex;
    }

}

void validateBranchCollection(std::vector<boost::shared_ptr<gridpack::component::DataCollection> > & branchCollectionVector,
    gridpack::parser::GOSSParser & parser)
{
    bool                     branchValid    = true;
    int                      branchIndex    = 0;

    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >::iterator dataCollection;
    for (dataCollection = branchCollectionVector.begin();
            dataCollection != branchCollectionVector.end(); ++dataCollection)
    {
        branchValid     = validateBranch(**dataCollection, parser, branchIndex);
        ++branchIndex;
//        parser.test_dumpDataColletion(**dataCollection);
    }
}

void validateCaseData(std::string caseID, int caseSBase)
{
    bool                     valid          = true;
    bool                     caseValid      = true;

    valid               = !strcmp("Case 1", caseID.c_str());
    caseValid           = valid;
    BOOST_CHECK_MESSAGE(valid,
        "CASE_ID actual                = --Case 1--  recovered is ++" << caseID << "++");

    valid               = caseSBase == 100;
    caseValid           = valid;
    BOOST_CHECK_MESSAGE(valid,
        "CASE_SBASE actual             = 100         recovered is " << caseSBase);

    BOOST_CHECK_MESSAGE(caseValid, "Case data was incorrect");
}



BOOST_AUTO_TEST_CASE(ArtificialData)
{
    bool ok = true;
    std::string   fileName("gridpack-test1.xml");
    std::cout << "TESTING LOAD OF " << fileName << std::endl;
    gridpack::parser::GOSSParser           parser;

    try {
        std::vector<boost::shared_ptr<gridpack::component::DataCollection> > busCollection;
        std::vector<boost::shared_ptr<gridpack::component::DataCollection> > branchCollection;
        parser.parse(fileName.c_str());

        std::string          caseID        = parser.getCaseId();
        int                  caseSBase     = parser.getCaseSbase();

        parser.copyDataCollection(busCollection, branchCollection);
//        parser.test_dumpDataColletionVector(branchCollection);
        validateBusCollection(busCollection, parser);
        validateBranchCollection(branchCollection, parser);
        validateCaseData(caseID, caseSBase);
    } catch (boost::property_tree::xml_parser_error & e) {
        e.what();
    } catch (boost::exception & e) {
        std::cout << "General exception\n\t";
        std::cout << std::endl;
    }

    std::cout << "END TESTING LOAD OF " << fileName << std::endl;


}

BOOST_AUTO_TEST_CASE(ArtificialData2)
{
    /*
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
    */
}

BOOST_AUTO_TEST_CASE(Greek_118)
{
/*
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
    */
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


