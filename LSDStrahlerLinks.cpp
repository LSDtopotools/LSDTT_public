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

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// create function: this makes a new StrahlerLinks object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
  
     //cout << "This stream order is: " << SO << " and the number of sources is: " << NThisOrderSources << endl;
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
      
      //cout <<"Order: " << SO-1 << " NS: " << SJunctions[SO-2].size() << " and NR: " << RJunctions[SO-2].size() << endl;
      
      // now sort and then loop through receivers to get the sources for the next
      // stream order.
      int NReceivers = int(thisOrderReceivers.size());
      sort(thisOrderReceivers.begin(), thisOrderReceivers.end());
      
      // reset the sources 
      thisOrderSources = emptyvec;
			
			// remove any receivers that have a stream order more than 1 greater than the order you are checking
			vector<int> NewReceivers;
			for (int r = 0; r < NReceivers; r++)
			{
				int ReceiverSO = JNetwork.get_StreamOrder_of_Junction(FlowInfo, thisOrderReceivers[r]);
				if (ReceiverSO == SO)
				{
          int same_SO = JNetwork.check_stream_order_of_upstream_nodes(thisOrderReceivers[r], FlowInfo);
          if (same_SO == 0)
          {
				    NewReceivers.push_back(thisOrderReceivers[r]);
          }
				}
			}
		
			NReceivers = int(NewReceivers.size());
			
      // get the starting receiver
      int LastReceiver;
      if (NReceivers > 0)
      {
				LastReceiver = NewReceivers[0];
        thisOrderSources.push_back(LastReceiver);
      }
          
      // go through the receivers removing duplicates
      if (NReceivers > 1)
      {
        for(int r = 1; r<NReceivers; r++)
        {
          // check to see if it is a new receiver
					if(NewReceivers[r] != LastReceiver)
					{
						LastReceiver = NewReceivers[r];
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

  // now get the nodes of the sources and recievers
  populate_NodeRowCol_vecvecs(JNetwork, FlowInfo);
                                
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to get a number of secondary data elements
// including the node indices of the downstream node in a link. 
// This is necessary since the downstream junction of a link is at a higher
// stream order and therefore not part of the individual link. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDStrahlerLinks::populate_NodeRowCol_vecvecs(LSDJunctionNetwork& JNetwork, 
                                                   LSDFlowInfo& FlowInfo)
{
  // temp vecvecs for holding the data
  vector< vector<int> > source_nodes;
  vector< vector<int> > source_rows;
  vector< vector<int> > source_cols;
  vector< vector<int> > receiver_nodes;
  vector< vector<int> > receiver_rows;
  vector< vector<int> > receiver_cols;

  vector<int> empty_vec;
  
  // these data elements hold local nodes, rows and columns
  vector<int> snodes,srows,scols,rnodes,rrows,rcols;
  int this_source_node;
  int this_source_row, this_source_col;
  int this_receiver_node;
  int this_receiver_row, this_receiver_col;
  int rj_node, last_receiver_node, last_receiver_row, last_receiver_col;

  // get the number of orders
  int NOrders = int(SourceJunctions.size());
  
  // loop through orders collecting data
  for (int order = 0; order<NOrders; order++)
  {
    // reset the local vectors
    snodes = empty_vec;
    srows = empty_vec;
    scols = empty_vec;
    rnodes = empty_vec;
    rrows = empty_vec;
    rcols = empty_vec;    
  
    // now loop through the starting and ending junctions
    int n_links_in_order = int(SourceJunctions[order].size());
    for(int link = 0; link<n_links_in_order; link++)
    {
      // get the source nodes
      this_source_node = JNetwork.get_Node_of_Junction(SourceJunctions[order][link]);
      FlowInfo.retrieve_current_row_and_col(this_source_node,this_source_row,this_source_col);  
    
      // now get the node of the receiver junction. This is one node downstream
      // of the terminating link node. 
      rj_node =  JNetwork.get_Node_of_Junction(ReceiverJunctions[order][link]);
      
      // look downstream until you hit the reciever node
      this_receiver_node = this_source_node;
      this_receiver_row = this_source_row;
      this_receiver_col = this_source_col;
      do
      {
        last_receiver_node = this_receiver_node;
        last_receiver_row = this_receiver_row;
        last_receiver_col = this_receiver_col;      
        FlowInfo.retrieve_receiver_information(last_receiver_node,this_receiver_node, 
                                             this_receiver_row, this_receiver_col);
      } while(this_receiver_node != rj_node);
      
      // populate the vectors
      snodes.push_back(this_source_node);
      srows.push_back(this_source_row);
      scols.push_back(this_source_col);
      rnodes.push_back(last_receiver_node);
      rrows.push_back(last_receiver_row);
      rcols.push_back(last_receiver_col);      
    }
    
    // insert the node row and column vectors into the vecvecs
    source_nodes.push_back(snodes);
    source_rows.push_back(srows);
    source_cols.push_back(scols);
    receiver_nodes.push_back(rnodes);
    receiver_rows.push_back(rrows);
    receiver_cols.push_back(rcols);
  }

  // copy the data members
  SourceNodes = source_nodes;
  SourceRows = source_rows;
  SourceCols = source_cols;
  ReceiverNodes = receiver_nodes;
  ReceiverRows = receiver_rows;
  ReceiverCols = receiver_cols;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calcualtes the drop of each link
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDStrahlerLinks::calculate_drops(LSDFlowInfo& FlowInfo, LSDRaster& topo_raster)
{
  vector< vector<float> > drops;
  vector<float> this_drop;
  vector<float> empty_vec;
  
  int srow,scol,rrow,rcol;
  float source_elev;
  float receiver_elev;
  
  // get the number of orders
  int NOrders = int(SourceJunctions.size());
  
  // loop through orders collecting data
  for (int order = 0; order<NOrders; order++)
  { 
    int n_links_in_order = int(SourceJunctions[order].size());
    
    // reset drop vector
    this_drop = empty_vec;
    
    for(int link = 0; link<n_links_in_order; link++)
    {
      srow = SourceRows[order][link];
      scol = SourceCols[order][link];
      rrow = ReceiverRows[order][link];
      rcol = ReceiverCols[order][link];
      
      // get the elevations of the source and receiver
      source_elev = topo_raster.get_data_element(srow,scol);  
      receiver_elev = topo_raster.get_data_element(rrow,rcol); 
      
      // add the drop to the drop vector
      this_drop.push_back(source_elev-receiver_elev);
      
    }
    
    // add the drop vector to the vecvec
    drops.push_back(this_drop);
  }
  
  DropData = drops;
    
}                                                   
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints the drops for assimilation into R or python
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDStrahlerLinks::print_drops(string data_directory, string threshold_string)
{
  int NOrders = int(DropData.size()); 
  string fname;
  string order_string;

  if (NOrders == 0)
  {
    cout << "You haven't calculated the drops yet! Not printing" << endl; 
  }
  else
  {
    for(int order = 0; order<NOrders; order++)
    {
      order_string = itoa(order);
      fname = data_directory+"Drops_Order_"+order_string+ "_Thresh_"+threshold_string;
      cout << "fname is: " << fname << endl;
      
      ofstream drops_out;
      drops_out.open(fname.c_str());
      
      int n_links_in_order = int(DropData[order].size());
    
      for(int link = 0; link<n_links_in_order; link++)
      {
        drops_out << DropData[order][link] << endl;
      }
      drops_out.close(); 
    
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns an LSDIndex raster of basins
// that do not receive drainage from nodes bordering the edge or nodes 
// bordering nodata
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDStrahlerLinks::get_no_edge_influence_mask(LSDFlowInfo& FI, 
                               LSDIndexRaster& Influence_Mask)
{
  // this array holds the data for the 
  Array2D<int> NoEdgeInfluence(NRows,NCols,0); 
  
  int IM_NDV = Influence_Mask.get_NoDataValue();
  
  // set the nodata points to nodata
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if(Influence_Mask.get_data_element(row,col) == IM_NDV)
      {
        NoEdgeInfluence[row][col] = IM_NDV;
      }
    }
  }
  
  // some parameters for holding data within the loop
  int this_row,this_col;
  vector<int> upslope_nodes;
  int this_us_row,this_us_col;
  
  // loop through all the drainage orders, 
  int NOrders = int(ReceiverNodes.size()); 

  // go in descending order, since you can skip some basins if they are 
  // in a non-affected channel
  for (int order = NOrders-1; order>=0; order--)
  {
    // get the number of links
    int NLinks = int(ReceiverNodes[order].size());
    for(int link = 0; link<NLinks; link++)
    {
      // get the location of the outlet
      this_row = ReceiverRows[order][link];
      this_col = ReceiverCols[order][link];
      
      // see if it is masked. If the Influence mask == 0, then it isn't masked
      // and we can include the basin in the valid data.
      if(Influence_Mask.get_data_element(this_row,this_col) == 0)
      {
        // check to see if this basin has already been tagged
        if(NoEdgeInfluence[this_row][this_col] != 1)
        {
          
          // it hasn't been tagged. Tag it. 
          upslope_nodes = FI.get_upslope_nodes(ReceiverNodes[order][link]); 
          
          // loop through the upslope nodes setting the NoEdgeInfuluence value to 1
          int Nupslope = int(upslope_nodes.size());
          for(int usn = 0; usn<Nupslope; usn++)
          {
            FI.retrieve_current_row_and_col(upslope_nodes[usn],this_us_row,
                                             this_us_col); 
            NoEdgeInfluence[this_us_row][this_us_col] = 1;
          }
        }
      }  
    } 
  }
  
  // now write the mask as an LSDIndexRaster
  LSDIndexRaster notInfluence_by_NDV(NRows,NCols,XMinimum,YMinimum,
                DataResolution,int(NoDataValue),NoEdgeInfluence,GeoReferencingStrings);
  return notInfluence_by_NDV;  

}                               
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function is a one stop function that returns an LSDRaster
// that has any pixel that has contributing pixels from the edge
// masked
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDStrahlerLinks::get_no_edge_influence_raster(LSDFlowInfo& FI,
                               LSDRaster& topography)
{
  // now look for the masked raster
  LSDIndexRaster mask = topography.find_cells_bordered_by_nodata();
  
  // get the influence mask
  LSDIndexRaster influence_mask = FI.find_cells_influenced_by_nodata(mask,topography);

  // now get the influence mask
  LSDIndexRaster NoEdgeInfluence = get_no_edge_influence_mask(FI,
                                                               influence_mask);

  // now get the masked raster
  int no_edge_influence_key = 0;
  LSDRaster masked_topography = topography.mask_to_nodata_with_mask_raster(NoEdgeInfluence, 
                                      no_edge_influence_key);

  // return this new raster
  return masked_topography;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the number of streams for each stream order
// FJC and MAH 17/03/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDStrahlerLinks::get_number_of_streams()
{
  vector<int> number_streams;
	int SO = 1;
  
  for (int i =0; i < int(SourceJunctions.size()); i++)
  {
    number_streams.push_back(int(SourceJunctions[i].size()));
    cout << "Stream order: " << SO << " N streams: " << SourceJunctions[i].size() << endl;
		SO++;
  }
  
  return number_streams;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#endif