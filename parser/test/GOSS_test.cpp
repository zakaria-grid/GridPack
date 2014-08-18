/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   PTI23_test.cpp
 * @author Kevin Glass
 * @date   2014-02-11 10:28:46 d3g096
 * 
 * @brief  Test PTI23_parser capability. Currently not implemented.
 * 
 * 
 */

#include <iostream>
#include <string>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <gridpack/parser/PTI23_parser.hpp>
#include <gridpack/configuration/configuration.hpp>
#include <gridpack/timer/coarse_timer.hpp>

#include "mpi.h"
#include <macdecls.h>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

// these values have to correspond to the actual values in the gridpack.xml file

#define EPSILON     0.0000001
#define TOLERANCE(x, y, eps) (x - EPSILON > y && x + EPSILON < y)

struct TestBranch
{
    int                      _BRANCH_FROMBUS;
    int                      _BRANCH_TOBUS;
    int                      elementIndex;
    int                      _BRANCH_INDEX;
    int                      mrid;
    double                   _BRANCH_FLOW_P;
    double                   _BRANCH_FLOW_Q;
    double                   _BRANCH_R;
    double                   _BRANCH_RATING_A;
    double                   _BRANCH_RATING_B;
    double                   _BRANCH_RATING_C;
    double                   _BRANCH_RATING;
    int                      _BRANCH_STATUS;
    double                   _BRANCH_X;
    double                   _BRANCH_B;
    double                   _BRANCH_SHUNT_ADMTTNC_BI;
    double                   _BRANCH_SHUNT_ADMTTNC_BJ;
    double                   _BRANCH_SHUNT_ADMTTNC_GI;
    double                   _BRANCH_SHUNT_ADMTTNC_GJ;
}
       /* bus data consists of:
         *    BUS_AREA
         *    BUS_BASEKV
         *    BUS_NAME
         *    BUS_NUMBER
         *    BUS_TYPE
         *    GENERATOR_BUSNUMBER
         *    GENERATOR_ID/
         *    GENERATOR_PG
         *    GENERATOR_QG
         *    GENERATOR_QMAX
         *    GENERATOR_QMIN
         *    GENERATOR_VS
         *    GENERATOR_IREG
         *    GENERATOR_MBASE
         *    GENERATOR_ZR
         *    GENERATOR_ZX
         *    GENERATOR_RT
         *    GENERATOR_XT
         *    GENERATOR_GTAP
         *    GENERATOR_STAT
         *    GENERATOR_RMPCT
         *    GENERATOR_PMAX
         *    GENERATOR_PMIN
         *    GENERATOR_OWNER
         *    BUS_LOAD_PL
         *    BUS_LOAD_QL
         *    LOAD_BUSNUMBER
         *    LOAD_STATUS
         *    LOAD_AREA
         *    LOAD_ZONE
         *    LOAD_PL
         *    LOAD_QL
         *    LOAD_IP
         *    LOAD_IQ
         *    LOAD_YP
         *    LOAD_YQ
         *    LOAD_OWNER
         *    BUS_OWNER
         *    powerGridId
         *    BUS_SHUNT_BL
         *    BUS_SHUNT_GL
         *    BUS_VOLTAGE_ANG
         *    BUS_VOLTAGE_MAG
         *    BUS_ZONE
         */
struct TestBus
{
    int                      _BUS_AREA;
    double                   _BUS_BASEKV;
    char *                   _BUS_NAME;
    int                      _BUS_NUMBER;
    int                      _BUS_TYPE;
    int                      _GENERATOR_BUSNUMBER;
    int                      _GENERATOR_ID;
    double                   _GENERATOR_PG;
    double                   _GENERATOR_QG;
    double                   _GENERATOR_QMAX;
    double                   _GENERATOR_QMIN;
    double                   _GENERATOR_VS;
    int                      _GENERATOR_IREG;
    double                   _GENERATOR_MBASE;
    double                   _GENERATOR_ZR;
    double                   _GENERATOR_ZX;
    double                   _GENERATOR_RT;
    double                   _GENERATOR_XT;
    double                   _GENERATOR_GTAP;
    int                      _GENERATOR_STAT;
    double                   _GENERATOR_RMPCT;
    double                   _GENERATOR_PMAX;
    double                   _GENERATOR_PMIN;
    char *                   _GENERATOR_OWNER;
    double                   _BUS__LOAD_PL;
    double                   _BUS__LOAD_QL;
    int                      _LOAD_BUSNUMBER;
    int                      _LOAD_STATUS;
    int                      _LOAD_AREA;
    int                      _LOAD_ZONE;
    double                   _LOAD_PL;
    double                   _LOAD_QL;
    double                   _LOAD_IP;
    double                   _LOAD_IQ;
    double                   _LOAD_YP;
    double                   _LOAD_YQ;
    char *                   _LOAD_OWNER;
    char *                   _BUS_OWNER;
    int                      powerGridId;
    double                   _BUS_SHUNT_BL;
    double                   _BUS_SHUNT_GL;
    double                   _BUS_VOLTAGE_ANG;
    double                   _BUS_VOLTAGE_MAG;
    int                      _BUS_ZONE;
};


BOOST_AUTO_TEST_SUITE(Parser)

void verifyNBranches(gridpack::network::BaseNetwork<TestBus, TestBranch> & network,
        int nBranches)
{
    bool                     correct        = true;
    if (network.numBranches() != nBranches) {
        printf("The number of branches should have been %d, but the network reported %d\n",network.numBranches(),nBranches);
        correct = false;
    }
    BOOST_CHECK_EQUAL(network.numBranches(), nBranches);
    BOOST_CHECK(correct);
}

void verifyNBuses(gridpack::network::BaseNetwork<TestBus, TestBranch> & network,
        int nBuses)
{
    bool                     correct        = true;
    if (network.numBuses() != nBuses) {
        printf("The number of buses should have been %d, but the network reported %d\n",network.numBuses(),nBuses);
        correct = false;
    }
    BOOST_CHECK_EQUAL(network.numBuses(), nBuses);
    BOOST_CHECK(correct);

}

void verifyBranchMaps(gridpack::network::BaseNetwork<TestBus, TestBranch> & network,
        int nMaps)
{
    /*
    bool                     correct        = true;
    BOOST_CHECK_EQUAL(network.numMap(), nMaps);
    BOOST_CHECK(correct);
    if (!correct) {
        printf("The number of branches should have been %d, but the network reported %d\n",
                    network.numMap(),nMaps);
        correct = false;
    }
    */
}

void verifyBranchData(gridpack::network::BaseNetwork<TestBus, TestBranch> & network, int branch, TestBus & testBranch)
{
    boost::shared_ptr<component::DataCollection>    branchCollection  = network.getBranch(branch);

    // get data collection object for "branch" in "network"
    BOOST_CHECK_EQUAL(testBranch._BRANCH_FROMBUS, branchCollection->getValue(BRANCH_FROMBUS));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: BRANCH_FROMBUS should have been %d, but was stored as %d\n", branch,
                testBranch._BRANCH_FROMBUS, branchCollection->getValue(BRANCH_FROMBUS));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_TOBUS, branchCollection->getValue(BRANCH_TOBUS));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d: BRANCH_TOBUS should have been %d, but was stored as %d\n", branch,
                testBranch._BRANCH_TOBUS, branchCollection->getValue(BRANCH_TOBUS));

    }

    BOOST_CHECK_EQUAL(testBranch.elementIndex, branchCollection->getValue(elementIndex));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: elementIndex should have been %d, but was stored as %d\n", branch,
                testBranch.elementIndex, branchCollection->getValue(elementIndex));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_INDEX, branchCollection->getValue(BRANCH_INDEX));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: BRANCH_INDEX should have been %d, but was stored as %d\n", branch,
                testBranch._BRANCH_INDEX, branchCollection->getValue(BRANCH_INDEX));

    }

    BOOST_CHECK_EQUAL(testBranch.mrid, branchCollection->getValue(mrid));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: mrid should have been %d, but was stored as %d\n", branch,
                testBranch.mrid, branchCollection->getValue(mrid));

    }


    BOOST_CHECK_EQUAL(testBranch._BRANCH_FLOW_P, branchCollection->getValue(BRANCH_FLOW_P));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: BRANCH_FLOW_P should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_FLOW_P, branchCollection->getValue(BRANCH_FLOW_P));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_FLOW_Q, branchCollection->getValue(BRANCH_FLOW_Q));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: BRANCH_FLOW_Q should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_FLOW_Q, branchCollection->getValue(BRANCH_FLOW_Q));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_R, branchCollection->getValue(BRANCH_R));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: BRANCH_R should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_R, branchCollection->getValue(BRANCH_R));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_RATING_A, branchCollection->getValue(BRANCH_RATING_A));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_RATING_A should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_RATING_A, branchCollection->getValue(BRANCH_RATING_A));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_RATING_B, branchCollection->getValue(BRANCH_RATING_B));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_RATING_B should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_RATING_B, branchCollection->getValue(BRANCH_RATING_B));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_RATING_C, branchCollection->getValue(BRANCH_RATING_C));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_RATING_C should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_RATING_C, branchCollection->getValue(BRANCH_RATING_C));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_RATING, branchCollection->getValue(BRANCH_RATING));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_RATING should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_RATING, branchCollection->getValue(BRANCH_RATING));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_STATUS, branchCollection->getValue(BRANCH_STATUS));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_STATUS should have been %d, but was stored as %d\n", branch,
                testBranch._BRANCH_STATUS, branchCollection->getValue(BRANCH_STATUS));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_X, branchCollection->getValue(BRANCH_X));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_X should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_X, branchCollection->getValue(BRANCH_X));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_B, branchCollection->getValue(BRANCH_B));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_B should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_B, branchCollection->getValue(BRANCH_B));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_SHUNT_ADMTTNC_BI, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_BI));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_SHUNT_ADMTTNC_BI should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_SHUNT_ADMTTNC_BI, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_BI));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_SHUNT_ADMTTNC_BJ, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_BJ));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_SHUNT_ADMTTNC_BJ should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_SHUNT_ADMTTNC_BJ, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_BJ));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_SHUNT_ADMTTNC_GI, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_GI));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Branch %d: _BRANCH_SHUNT_ADMTTNC_GI should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_SHUNT_ADMTTNC_GI, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_GI));

    }

    BOOST_CHECK_EQUAL(testBranch._BRANCH_SHUNT_ADMTTNC_GJ, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_GJ));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d: _BRANCH_SHUNT_ADMTTNC_GJ should have been %f, but was stored as %f\n", branch,
                testBranch._BRANCH_SHUNT_ADMTTNC_GJ, branchCollection->getValue(BRANCH_SHUNT_ADMTTNC_GJ));

    }
}

void verifyBusData(gridpack::network::BaseNetwork<TestBus, TestBranch> & network, int bus, TestBus & testBus)
{
    bool                      correct       = true;
    boost::shared_ptr<component::DataCollection>  busCollection  = network.getBus(bus);
    
    // compare 
    BOOST_CHECK_EQUAL(testBus._BUS_AREA, busCollection->getValue(BUS_AREA));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_AREA should have been %d, but was stored as %d\n", bus,
                testBus._BUS_AREA, busCollection->getValue(BUS_AREA));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_BASEKV, busCollection->getValue(BUS_BASEKV));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_BASEKV should have been %f, but was stored as %f\n", bus,
                testBus._BUS_BASEKV, busCollection->getValue(BUS_BASEKV));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_NAME, busCollection->getValue(BUS_NAME));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_BUS_NAME should have been %s, but was stored as %c\n", bus,
                testBus._BUS_NAME, busCollection->getValue(BUS_NAME));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_NUMBER, busCollection->getValue(BUS_NUMBER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_NUMBER should have been %d, but was stored as %d\n", bus,
                testBus._BUS_NUMBER, busCollection->getValue(BUS_NUMBER));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_TYPE, busCollection->getValue(BUS_TYPE));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_BUS_TYPE should have been %d, but was stored as %d\n", bus,
                testBus._BUS_TYPE, busCollection->getValue(BUS_TYPE));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_BUSNUMBER, busCollection->getValue(GENERATOR_BUSNUMBER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_BUSNUMBER should have been %d, but was stored as %d\n", bus,
                testBus._GENERATOR_BUSNUMBER, busCollection->getValue(GENERATOR_BUSNUMBER));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_ID, busCollection->getValue(GENERATOR_ID));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_GENERATOR_ID should have been %d, but was stored as %d\n", bus,
                testBus._GENERATOR_ID, busCollection->getValue(GENERATOR_ID));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_PG, busCollection->getValue(GENERATOR_PG));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_PG should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_PG, busCollection->getValue(GENERATOR_PG));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_QG, busCollection->getValue(GENERATOR_QG));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_QG should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_QG, busCollection->getValue(GENERATOR_QG));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_QMAX, busCollection->getValue(GENERATOR_QMAX));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_QMAX should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_QMAX, busCollection->getValue(GENERATOR_QMAX));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_QMIN, busCollection->getValue(GENERATOR_QMIN));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_QMIN should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_QMIN, busCollection->getValue(GENERATOR_QMIN));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_VS, busCollection->getValue(GENERATOR_VS));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_VS should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_VS, busCollection->getValue(GENERATOR_VS));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_IREG, busCollection->getValue(GENERATOR_IREG));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_IREG should have been %d, but was stored as %d\n", bus,
                testBus._GENERATOR_IREG, busCollection->getValue(GENERATOR_IREG));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_MBASE, busCollection->getValue(GENERATOR_MBASE));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_MBASE should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_MBASE, busCollection->getValue(GENERATOR_MBASE));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_ZR, busCollection->getValue(GENERATOR_ZR));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_ZR should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_ZR, busCollection->getValue(GENERATOR_ZR));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_ZX, busCollection->getValue(GENERATOR_ZX));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_ZX should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_ZX, busCollection->getValue(GENERATOR_ZX));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_RT, busCollection->getValue(GENERATOR_RT));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_RT should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_RT, busCollection->getValue(GENERATOR_RT));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_XT, busCollection->getValue(GENERATOR_XT));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_XT should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_XT, busCollection->getValue(GENERATOR_XT));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_GTAP, busCollection->getValue(GENERATOR_GTAP));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_GTAP should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_GTAP, busCollection->getValue(GENERATOR_GTAP));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_STAT, busCollection->getValue(GENERATOR_STAT));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_STAT should have been %d, but was stored as %d\n", bus,
                testBus._GENERATOR_STAT, busCollection->getValue(GENERATOR_STAT));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_RMPCT, busCollection->getValue(GENERATOR_RMPCT));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_RMPCT should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_RMPCT, busCollection->getValue(GENERATOR_RMPCT));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_PMAX, busCollection->getValue(GENERATOR_PMAX));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_PMAX should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_PMAX, busCollection->getValue(GENERATOR_PMAX));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_PMIN, busCollection->getValue(GENERATOR_PMIN));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_PMIN should have been %f, but was stored as %f\n", bus,
                testBus._GENERATOR_PMIN, busCollection->getValue(GENERATOR_PMIN));

    }

    BOOST_CHECK_EQUAL(testBus._GENERATOR_OWNER, busCollection->getValue(GENERATOR_OWNER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _GENERATOR_OWNER should have been %s, but was stored as %c\n", bus,
                testBus._GENERATOR_OWNER, busCollection->getValue(GENERATOR_OWNER));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_LOAD_PL, busCollection->getValue(BUS_LOAD_PL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_LOAD_PL should have been %f, but was stored as %f\n", bus,
                testBus._BUS_LOAD_PL, busCollection->getValue(BUS_LOAD_PL));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_LOAD_QL, busCollection->getValue(BUS_LOAD_QL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_LOAD_QL should have been %f, but was stored as %f\n", bus,
                testBus._BUS_LOAD_QL, busCollection->getValue(BUS_LOAD_QL));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_BUSNUMBER, busCollection->getValue(LOAD_BUSNUMBER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_LOAD_BUSNUMBER should have been %d, but was stored as %d\n", bus,
                testBus._LOAD_BUSNUMBER, busCollection->getValue(LOAD_BUSNUMBER));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_STATUS, busCollection->getValue(LOAD_STATUS));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_STATUS should have been %d, but was stored as %d\n", bus,
                testBus._LOAD_STATUS, busCollection->getValue(LOAD_STATUS));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_AREA, busCollection->getValue(LOAD_AREA));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_LOAD_AREA should have been %d, but was stored as %d\n", bus,
                testBus._LOAD_AREA, busCollection->getValue(LOAD_AREA));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_ZONE, busCollection->getValue(LOAD_ZONE));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_ZONE should have been %d, but was stored as %d\n", bus,
                testBus._LOAD_ZONE, busCollection->getValue(LOAD_ZONE));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_PL, busCollection->getValue(LOAD_PL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_PL should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_PL, busCollection->getValue(LOAD_PL));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_QL, busCollection->getValue(LOAD_QL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_QL should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_QL, busCollection->getValue(LOAD_QL));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_IP, busCollection->getValue(LOAD_IP));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_IP should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_IP, busCollection->getValue(LOAD_IP));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_IQ, busCollection->getValue(LOAD_IQ));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_IQ should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_IQ, busCollection->getValue(LOAD_IQ));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_YP, busCollection->getValue(LOAD_YP));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d :_LOAD_YP should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_YP, busCollection->getValue(LOAD_YP));

    }

    BOOST_CHECK_EQUAL(testBus._LOAD_YQ, busCollection->getValue(LOAD_YQ));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_YQ should have been %f, but was stored as %f\n", bus,
                testBus._LOAD_YQ, busCollection->getValue(LOAD_YQ));

    }


    BOOST_CHECK_EQUAL(testBus._LOAD_OWNER, busCollection->getValue(LOAD_OWNER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _LOAD_OWNER should have been %s, but was stored as %c\n", bus,
                testBus._LOAD_OWNER, busCollection->getValue(LOAD_OWNER));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_OWNER, busCollection->getValue(BUS_OWNER));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_OWNER should have been %s, but was stored as %c\n", bus,
                testBus._BUS_OWNER, busCollection->getValue(BUS_OWNER));

    }

    BOOST_CHECK_EQUAL(testBus.powerGridId, busCollection->getValue(powerGridId));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : powerGridId should have been %d, but was stored as %d\n", bus,
                testBus.powerGridId, busCollection->getValue(powerGridId));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_SHUNT_BL, busCollection->getValue(BUS_SHUNT_BL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_SHUNT_BL should have been %f, but was stored as %f\n", bus,
                testBus._BUS_SHUNT_BL, busCollection->getValue(BUS_SHUNT_BL));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_SHUNT_GL, busCollection->getValue(BUS_SHUNT_GL));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_SHUNT_GL should have been %f, but was stored as %f\n", bus,
                testBus._BUS_SHUNT_GL, busCollection->getValue(BUS_SHUNT_GL));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_VOLTAGE_ANG, busCollection->getValue(BUS_VOLTAGE_ANG));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_VOLTAGE_ANG should have been %f, but was stored as %f\n", bus,
                testBus._BUS_VOLTAGE_ANG, busCollection->getValue(BUS_VOLTAGE_ANG));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_VOLTAGE_MAG, busCollection->getValue(BUS_VOLTAGE_MAG));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_VOLTAGE_MAG should have been %f, but was stored as %f\n", bus,
                testBus._BUS_VOLTAGE_MAG, busCollection->getValue(BUS_VOLTAGE_MAG));

    }

    BOOST_CHECK_EQUAL(testBus._BUS_ZONE, busCollection->getValue(BUS_ZONE));
    BOOST_CHECK(correct);
    if (!correct) {
        printf("Bus %d : _BUS_ZONE should have been %d, but was stored as %d\n", bus,
                testBus._BUS_ZONE, busCollection->getValue(BUS_ZONE));

    }
}

/*
 * Load a network with the data from gridpack-test.xml
 * Verify nBranches and nBuses
 * Verify transmission elements between branches
 * Use the network to retreive data from
 */
BOOST_AUTO_TEST_CASE( DataCompletenessTest )
{

    int ierr;
    int me;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    int nprocs;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Create network
    gridpack::parallel::Communicator world;
    gridpack::network::BaseNetwork<TestBus, TestBranch> network(world);
    gridpack::parser::PTI23_parser<gridpack::network::BaseNetwork<TestBus, TestBranch> parser(network);

    // only let proc 0 read
    if (me == 0) {
        bool ok = true;
        int             n                   = NBUSES;

        std::string   fileName("gridpack-test1.xml");
        parser.parse(fileName);

        /*
         * For gridpack-test.xml
         * The branches are:
         *  Branch      x->y
         *  1           1->2
         *  2           3->4
         *  3           5->6
         *  4           7->8
         *  5           9->10
         *  6           2->1
         *  7           1->2
         *  8           7->8
         *  9           2->10
         *  10          10->2
         *
         * The maps are:
         * Connection       Number          Switched
         *                  in
         *                  Connection
         *  1->2            3               1
         *  3->4            1               0
         *  5->6            1               0
         *  7->8            2               0
         *  9->10           1               0
         *  10->2           2               1
         */
        const int       NBRANCHES           = 10;
        const int       NBUSES              = 10;
        const int       NMAPS               = 6;

        verifyNBranches(network, 10);
        verifyNBuses(network, 10);
        verifyBranchMaps(network, 10);


        TestBus  testBus;

        for (int i = 0; i < NBUSES; i++) {
            testBus =
            {
                    1, i + .02, i + .03, i + .04, i + .05, i + .06, i, i + .07, i + .08, i + .09, i + .10, i + .11,
                    i + .12, i, i + .13, i + .14, i + .15, i, i + .16, i + .17, i, 1, i, 1, i + .18, i + .19,
                    i + .20, i + .21, i + .22, i + .23, i, 1, i, i + .24, i + .25, i + .26, i + .27
            };
            verifyBusData(network, 1,);
        }

        TestBus  testBus;

        TestBranch           testBranch;
        testBranch =
        {
                1, 2, 1, 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1, 1.08, 1.09,
                1.1, 1.2, 1.3, 1.4
        }
        verifyBranchData(network, 1, testBranch);
        testBranch =
        {
                3, 4, 1, 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1, 1.08, 1.09,
                2.1, 2.2, 2.3, 2.4
        }
        verifyBranchData(network, 2, testBranch);
        testBranch =
        {
                5, 6, 3, 3, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3, 3.08, 3.09,
                3.1, 3.2, 3.3, 3.4
        }
        verifyBranchData(network, 3, testBranch);
        testBranch =
        {
                7, 8, 4, 4, 4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4, 4.08, 4.09,
                4.1, 4.2, 4.3, 4.4
        }
        verifyBranchData(network, 4, testBranch);
        testBranch =
        {
                9, 10, 5, 5, 5.01, 5.02, 5.03, 5.04, 5.05, 5.06, 5.07, 5, 5.08, 5.09,
                5.1, 5.2, 5.3, 5.4
        }
        verifyBranchData(network, 5, testBranch);
        testBranch =
        {
                2, 1, 6, 6, 6.01, 6.02, 6.03, 6.04, 6.05, 6.06, 6.07, 6, 6.08, 6.09,
                6.1, 6.2, 6.3, 6.4
        }
        verifyBranchData(network, 6, testBranch);
        testBranch =
        {
                1, 2, 7, 7, 7.01, 7.02, 7.03, 7.04, 7.05, 7.06, 7.07, 7, 7.08, 7.09,
                7.1, 7.2, 7.3, 7.4
        }
        verifyBranchData(network, 7, testBranch);
        testBranch =
        {
                7, 8, 8, 8, 8.01, 8.02, 8.03, 8.04, 8.05, 8.06, 8.07, 8, 8.08, 8.09,
                8.1, 8.2, 8.3, 8.4
        }
        verifyBranchData(network, 8, testBranch);
        testBranch =
        {
                10, 2, 9, 9, 9.01, 9.02, 9.03, 9.04, 9.05, 9.06, 9.07, 9, 9.08, 9.09,
                9.1, 9.2, 9.3, 9.4
        }
        verifyBranchData(network, 9, testBranch);
        {
                2, 10, 10, 10, 10.01, 10.02, 10.03, 10.04, 10.05, 10.06, 10.07, 10, 10.08, 10.09,
                10.1, 10.2, 10.3, 10.4
        }
        verifyBranchData(network, 10, testBranch);
}

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
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  int me = world.rank();
  if (me == 0) {
    printf("Testing Serial Input\n");
  }
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  return result;
}


