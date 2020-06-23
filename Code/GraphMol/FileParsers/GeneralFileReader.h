//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef GENERAL_FILE_READER_H
#define GENERAL_FILE_READER_H
#include <string>
#include <iostream>
#include <vector>
#include "MolSupplier.h"
#include <RDGeneral/BadFileException.h>
#include <RDStreams/streams.h>

namespace RDKit{
	struct SupplierOptions{
		bool takeOwnership = true;
		bool sanitize = true;
		bool removeHs = true;
		bool strictParsing = true;

		std::string delimiter = "\t";
		int smilesColumn = 0;
		int nameColumn = 1;
		bool titleLine = true;

		std::string nameRecord = "";
		int confId2D = -1;
		int confId3D = 0;
	};

  class GeneralFileReader{
	  public:
			std::string d_path;
		  std::string d_fileFormat, d_compressionFormat;
		  std::vector<std::string> d_fileFormats{ "sdf", "mae", "smi", "csv", "txt", "tsv", "tdt"};
		  std::vector<std::string> d_compressionFormats{ "gz" };
		struct SupplierOptions d_opt;	//! d_options for the Mol Supplier 
		public:
			GeneralFileReader(const std::string& path);
		  GeneralFileReader(const std::string& path, const struct SupplierOptions options);
			//! Function to get the file name from the path
			std::string getFileName();
		  //! Function to check the validity of the file and compression format
		  bool valid();
		  //! Function to set the file and compression format
		  void determineFormat();
		  //! Get MolSupplier Object based on the file and compression format
 			/*!
    		<b>Note:</b> the caller is responsible for <tt>delete</tt>ing the result
  		*/
		  RDKit::MolSupplier* getSupplier();
  };

	GeneralFileReader::GeneralFileReader(const std::string& path){
		d_path = path;
		d_fileFormat = "";
		d_compressionFormat = "";
	}

	GeneralFileReader::GeneralFileReader(const std::string& path, const struct SupplierOptions options){
		d_path = path;
		d_fileFormat = "";
		d_compressionFormat = "";
		d_opt = options;
	}

	bool GeneralFileReader::valid(){
    //! Case 1: Unconventional format types
		if(d_fileFormat.compare("maegz") == 0){
			d_fileFormat = "mae";
			d_compressionFormat = "gz";
			return true;
		}

	
		//! Case 2: Either the filename has a file format or compression format or both
		bool flag_fileFormat = std::find(d_fileFormats.begin(), d_fileFormats.end(), d_fileFormat) != d_fileFormats.end();
		
		if (!d_compressionFormat.empty()){
			bool flag_compressionFormat = std::find(d_compressionFormats.begin(), d_compressionFormats.end(), d_compressionFormat) 
			                              != d_compressionFormats.end();
			                              
			//! if the compression type is not valid then 
			if(!flag_compressionFormat){
				d_fileFormat = d_compressionFormat;
				d_compressionFormat = "";
				flag_fileFormat = std::find(d_fileFormats.begin(), d_fileFormats.end(), d_fileFormat) != d_fileFormats.end();
			}	
    }
		return flag_fileFormat;
	}


	std::string RDKit::GeneralFileReader::getFileName(){
		
		char delimiter = '/';
		std::string fname  = "";
		auto slash = d_path.rfind(delimiter, d_path.length());
		if(slash != std::string::npos){
			fname += d_path.substr(slash + 1, d_path.length() - slash);
		} 	
		return fname;
	}

	void RDKit::GeneralFileReader::determineFormat(){
	
		std::string fileName = getFileName();
		int dots = std::count(fileName.begin(), fileName.end(), '.');

		if (dots == 0) throw std::invalid_argument("Recieved Invalid File Format, no extension or compression");

		else if (dots == 1){
			//! there is a file format but no compression format
			int pos = fileName.find(".");
			d_fileFormat = fileName.substr(pos + 1);
			if (!valid()) throw std::invalid_argument("Recieved Invalid File Format");
		}
		else{
			//! there is a file and compression format
			int n = fileName.length();
			int p1 = fileName.rfind(".");
			int p2 = fileName.rfind(".", p1 - 1);
			d_fileFormat = fileName.substr(p2 + 1, (p1 - p2) - 1);
			//! possible compression format
			d_compressionFormat = fileName.substr(p1 + 1, (n - p1) + 1);
			if (!valid()) throw std::invalid_argument("Recieved Invalid File or Compression Format");
		}
	}

	RDKit::MolSupplier *RDKit::GeneralFileReader::getSupplier(){

    determineFormat();
		if (d_compressionFormat.empty()){
			std::ifstream *strm = new std::ifstream(d_path.c_str());
			if (d_fileFormat.compare("sdf") == 0){
				ForwardSDMolSupplier *sdsup = new ForwardSDMolSupplier(strm, d_opt.takeOwnership, d_opt.sanitize, d_opt.removeHs, d_opt.strictParsing);
				return sdsup;
			}

			else if (d_fileFormat == "smi"){
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "csv"){
				d_opt.delimiter = ",";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "txt"){
				d_opt.delimiter = "\t";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "tsv"){
				d_opt.delimiter = "\t";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}
 
			else if (d_fileFormat == "mae"){
				MaeMolSupplier *maesup = new MaeMolSupplier(strm, d_opt.takeOwnership, d_opt.sanitize, d_opt.removeHs);
				return maesup;
			}

			else if (d_fileFormat == "tdt"){
				TDTMolSupplier *tdtsup = new TDTMolSupplier(strm, d_opt.takeOwnership, d_opt.nameRecord, d_opt.confId2D, d_opt.confId3D, d_opt.sanitize);
				return tdtsup;
			}
		} 
		else{
			auto *strm = new gzstream(d_path);

			if (d_fileFormat == "sdf"){
				ForwardSDMolSupplier *sdsup = new ForwardSDMolSupplier(strm, d_opt.takeOwnership, d_opt.sanitize, d_opt.removeHs, d_opt.strictParsing);
				return sdsup;
			}

			else if (d_fileFormat == "smi"){
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "csv"){
				d_opt.delimiter = ",";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "txt"){
				d_opt.delimiter = "\t";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}

			else if (d_fileFormat == "tsv"){
				d_opt.delimiter = "\t";
				SmilesMolSupplier *smsup = new SmilesMolSupplier(strm, d_opt.takeOwnership, d_opt.delimiter, d_opt.smilesColumn, d_opt.nameColumn, d_opt.titleLine, d_opt.sanitize);
				return smsup;
			}
 
			else if (d_fileFormat == "mae"){
				MaeMolSupplier *maesup = new MaeMolSupplier(strm, d_opt.takeOwnership, d_opt.sanitize, d_opt.removeHs);
				return maesup;
			}

			else if (d_fileFormat == "tdt"){
				TDTMolSupplier *tdtsup = new TDTMolSupplier(strm, d_opt.takeOwnership, d_opt.nameRecord, d_opt.confId2D, d_opt.confId3D, d_opt.sanitize);
				return tdtsup;
			}
		}
	}
}
#endif
