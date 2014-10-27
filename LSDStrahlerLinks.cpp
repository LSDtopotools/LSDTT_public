//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStrahlerLinks
// Land Surface Dynamics StrahlerLinks
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  that keeps track of Strahler ordered stream links and computes
//  various statistics on them
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2014 Simon M. Mudd 2014
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStrahlerLinks.cpp
// LSDStrahlerLinks object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.1.0		26/10/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDStrahlerLinks_CPP
#define LSDStrahlerLinks_CPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
//#include "LSDRaster.hpp"
//#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
//#include "LSDIndexChannel.hpp"
//#include "LSDStatsTools.hpp"
#include "LSDStrahlerLinks.hpp"
using namespace std;
using namespace TNT;


void LSDStrahlerLinks::create()
{
  cout << "You need to designate a Junction Network and a Flow info to initialise"
       << " a LSDStrahlerLinks object" << endl;
}

// create function: this makes a new StrahlerLinks object
void LSDStrahlerLinks::create(LSDJunctionNetwork& JNetwork, LSDFlowInfo& FlowInfo)
{
  // these data memebers will be replaced later
  vector< vector<int> >  SJunctions;
  vector< vector<int> >  RJunctions;
  
  // get georeferencing information
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();
  
  // now get all the sources and loop through them to get the 1st order basins
  vector<int> thisOrderNISources = JNetwork.get_SourcesVector();
  int NThisOrderSources = int(thisOrderNISources.size());  
  vector<int> thisOrderSources(NThisOrderSources);
  
  // convert these sources to junctions
  for(int i = 0; i<NThisOrderSources; i++)
  {
    thisOrderSources[i] = JNetwork.get_Junction_of_Node(thisOrderNISources[i],FlowInfo);
  }
  //cout << "Strahler, LINE 128 N_Sources: " <<  NThisOrderSources << endl;
  
  
  int ds_SO_link;
  vector<int> emptyvec;
  vector<int> thinnedSources;
  vector<int> thisOrderReceivers;
    

  
  int SO = 1;
  // go into a loop that keeps making source and receiver vectors until you 
  // run out of sources
  while(NThisOrderSources>0)
  {
  
     cout << "This stream order is: " << SO << " and the number of sources is: " 
          << NThisOrderSources << endl;
     SO++;
  
    // reset the sources and receivers
    thinnedSources = emptyvec;
    thisOrderReceivers = emptyvec;
    
    // loop through this order
    for(int s = 0; s<NThisOrderSources; s++)
    {
      // find the downstream link  
      ds_SO_link = JNetwork.get_Next_StreamOrder_Junction(thisOrderSources[s]);
      
      // discard sources that end in a baselevel node
      if(ds_SO_link != NoDataValue)
      {
        thinnedSources.push_back(thisOrderSources[s]);
        thisOrderReceivers.push_back(ds_SO_link);
      }     
    }
    
    if (int(thinnedSources.size()) > 0)
    {
      // now add the sources and receivers to the data memebers
      SJunctions.push_back(thinnedSources);
      RJunctions.push_back(thisOrderReceivers);
      
      cout <<"Order: " << SO-1 << " NS: " << SJunctions[SO-2].size() << " and NR: " << RJunctions[SO-2].size() << endl;
      
      // now sort and then loop through receivers to get the sources for the next
      // stream order.
      int NReceivers = int(thisOrderReceivers.size());
      sort(thisOrderReceivers.begin(), thisOrderReceivers.end());
      
      // reset the sources 
      thisOrderSources = emptyvec;
      
      // get the starting receiver
      int LastReceiver;
      if (NReceivers > 0)
      {
        LastReceiver = thisOrderReceivers[0];
        thisOrderSources.push_back(LastReceiver);
      }
          
      // go through the receivers removing duplicates
      if (NReceivers > 1)
      {
        for(int r = 1; r<NReceivers; r++)
        {
          // check to see if it is a new receiver
          if(thisOrderReceivers[r] != LastReceiver)
          {
            LastReceiver = thisOrderReceivers[r];
            thisOrderSources.push_back(LastReceiver);
          }                  
        }
      }
      
      NThisOrderSources = int(thisOrderSources.size());
    }
    else
    {
      NThisOrderSources = 0;  
    }    
  }
  
    
  // update data members
  SourceJunctions = SJunctions;
  ReceiverJunctions = RJunctions;
  
  // now get the number of stream orders
  int NOrders = int(SourceJunctions.size());
  int t_sources = 0;
  // make sure all the vectors are the right size
  for(int o = 0; o<NOrders; o++)
  {
    t_sources += SourceJunctions[o].size();
    cout << "order: " << o+1 << " nsources: " << SourceJunctions[o].size()
         << " and n receivers: " << ReceiverJunctions[o].size() << endl;
  }
  cout << "Total sources: " << t_sources << endl;
  
}

#endif