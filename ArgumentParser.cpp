#include"ArgumentParser.h"
#include<algorithm>
#include<cstdlib>
#include<sstream>

#include "common.h"

#define FF first
#define SS second
#define Sof(x) (long)x.size()

vector <double> tokenizeD(const string &input,const string &space = ","){//{{{
   vector <double> ret;
   long pos=0,f=0,n=input.size();
   while((pos<n)&&(f<n)&&(f>=0)){
      f=input.find(space,pos);
      if(f==pos)pos++;
      else{
         if((f <n)&&(f>=0)){
            ret.push_back(atof(input.substr(pos,f-pos).c_str()));
            pos=f+1;
         }
      }
   }
   if(pos<n)ret.push_back(atof(input.substr(pos,n-pos).c_str()));
   return ret;
} //}}}


// GET {{{
vector<string>& ArgumentParser::args(){
   return arguments;
}
string ArgumentParser::getS(string name){
   if(!existsOption(name))error("ArgumentParser: argument name %s unknown.\n",name.c_str());
   if(mapS.find(name)!=mapS.end())
      return mapS[name];
   return "";
}
long ArgumentParser::getL(string name){
   if(!existsOption(name))error("ArgumentParser: argument name %s unknown.\n",(name).c_str());
   if(mapL.find(name)!=mapL.end())
      return mapL[name];
   return -1;
}
double ArgumentParser::getD(string name){
   if(!existsOption(name))error("ArgumentParser: argument name %s unknown.\n",(name).c_str());
   if(mapD.find(name)!=mapD.end())
      return mapD[name];
   return -1;
}
bool ArgumentParser::flag(string name){
   if(!existsOption(name))error("ArgumentParser: argument name %s unknown.\n",(name).c_str());
   return isSet(name);
}//}}}
vector<double> ArgumentParser::getTokenizedS2D(string name){
   if(!existsOption(name))error("ArgumentParser: argument name %s unknown.\n",name.c_str());
   if(mapS.find(name)!=mapS.end())
      return tokenizeD(mapS[name]);
   return vector<double>();
}
bool ArgumentParser::parse(int argc,char * argv[]){//{{{
//   for(long i=0;i<argc;i++)message("_%s_\n",(args[i]).c_str());
   // add verbose if  possible {{{
   if(! (existsName("v")||existsName("verbose")||existsOption("verbose")))
      addOptionB("v","verbose","verbose",0,"Verbose output.");
   //if(! (existsName("V")||existsName("veryVerbose")||existsOption("veryVerbose")))
   //   addOptionB("V","veryVerbose","veryVerbose",0,"Very verbose output.");
   // }}}
   programName=(string)argv[0];
   string val,opt;
   for(long i = 1; i<argc;i++){
      val=(string)argv[i];
      if(val[0]!='-'){
         arguments.push_back(val);
         continue;
      }
      if(Sof(val)==2){
         opt=val.substr(1,1);
         val="";
      }else{
         if(val.find("=")!=string::npos){
            opt=val.substr(2,val.find("=")-2);
            val=val.substr(val.find("=")+1);
         }else{
            opt=val.substr(2);
            val="";
         }
      }
      if((opt=="help")||(opt=="h")){
         usage();
         return false;
      }
      if(names.find(opt)==names.end()){
         error("Unknown option '%s'.\n",argv[i]);
         return false;
      }
      if(validOptions[names[opt]].type!=OTBOOL){
         if(val==""){
            i++;
            if(i==argc)break;
            val = argv[i];
         }
         switch(validOptions[names[opt]].type){
            case OTSTRING:
               mapS[names[opt]]=val;
               break;
            case OTLONG:
               mapL[names[opt]]=atoi(val.c_str());
               break;
            case OTDOUBLE:
               mapD[names[opt]]=atof(val.c_str());
               break;
            case OTBOOL:;
         }
      }else{
         mapB[names[opt]]=!mapB[names[opt]];
      }
   }
   //writeAll();
   if(Sof(arguments)<minimumArguments){
      error("Need at least %ld arguments.\n\n",minimumArguments);
      usage();
      return false;
   }
   for(long i = 0;i<Sof(compulsory);i++){
      if(! isSet(compulsory[i])){
         error("Missing option \"%s\"\n",(compulsory[i]).c_str());
         usage();
         return false;
      }
   }
   // set public variable verbose 
   verbose = flag("verbose")||(existsOption("veryVerbose")&&flag("veryVerbose"));
   return true;
}//}}}
void ArgumentParser::writeAll(){//{{{
   message("arguments: ");
   for(long i=0;i<Sof(arguments);i++)
      message("%s ",(arguments[i]).c_str());
   message("\n");
   for(map<string,string>::iterator it=mapS.begin();it!=mapS.end();it++){
      message("OPT:%s VAL:%s\n",(it->FF).c_str(),(it->SS).c_str());
   }
   for(map<string,long>::iterator it=mapL.begin();it!=mapL.end();it++){
      message("OPT:%s VAL:%ld\n",(it->FF).c_str(),it->SS);
   }
   for(map<string,double>::iterator it=mapD.begin();it!=mapD.end();it++){
      message("OPT:%s VAL:%lf\n",(it->FF).c_str(),(it->SS));
   }
   for(map<string,bool>::iterator it=mapB.begin();it!=mapB.end();it++){
      message("OPT:%s VAL:%d\n",(it->FF).c_str(),(it->SS));
   }
}//}}}
void ArgumentParser::addOptionL(string shortName,string longName, string name, bool comp, string description, long defValue){//{{{
   Option newOpt;
   if(existsOption(name)){
      error("ArgumentParser: Option \"%s\"\n",(name).c_str());
      return;
   }
   if(defValue!=-47){
      appendDescription<long>(description,defValue);
      mapL[name]=defValue;
   }
   newOpt.type=OTLONG;
   newOpt.shortName=shortName;
   newOpt.longName=longName;
   newOpt.description=description;
   validOptions[name]=newOpt;
   if(shortName!="")names[shortName]=name;
   if(longName!="")names[longName]=name;
   if(comp)compulsory.push_back(name);
}//}}}
void ArgumentParser::addOptionD(string shortName,string longName, string name, bool comp, string description, double defValue){//{{{
   Option newOpt;
   if(existsOption(name)){
      error("ArgumentParser: Option \"%s\"\n",(name).c_str());
      return;
   }
   if(defValue!=-47.47){
      appendDescription<double>(description,defValue);
      mapD[name]=defValue;
   }
   newOpt.type=OTDOUBLE;
   newOpt.shortName=shortName;
   newOpt.longName=longName;
   newOpt.description=description;
   validOptions[name]=newOpt;
   if(shortName!="")names[shortName]=name;
   if(longName!="")names[longName]=name;
   if(comp)compulsory.push_back(name);
}//}}}
void ArgumentParser::addOptionB(string shortName,string longName, string name, bool comp, string description, bool defValue){//{{{
   Option newOpt;
   if(existsOption(name)){
      error("ArgumentParser: Option \"%s\"\n",(name).c_str());
      return;
   }
   mapB[name]=defValue;
   if(defValue) description +=" (default: On)";
   else description+=" (default: Off)";
   newOpt.type=OTBOOL;
   newOpt.shortName=shortName;
   newOpt.longName=longName;
   newOpt.description=description;
   validOptions[name]=newOpt;
   if(shortName!="")names[shortName]=name;
   if(longName!="")names[longName]=name;
   if(comp)compulsory.push_back(name);
}//}}}
void ArgumentParser::addOptionS(string shortName,string longName, string name, bool comp, string description, string defValue){//{{{
   Option newOpt;
   if(existsOption(name)){
      error("ArgumentParser: Option \"%s\"\n",(name).c_str());
      return;
   }
   if(defValue!="noDefault"){
      appendDescription<string>(description,defValue);
      mapS[name]=defValue;
   }
   newOpt.type=OTSTRING;
   newOpt.shortName=shortName;
   newOpt.longName=longName;
   newOpt.description=description;
   validOptions[name]=newOpt;
   if(shortName!="")names[shortName]=name;
   if(longName!="")names[longName]=name;
   if(comp)compulsory.push_back(name);
}//}}}
template <typename valueType>
void ArgumentParser::appendDescription(string &desc,valueType defValue){//{{{
   stringstream descStream;
   descStream<<desc<<" (default: "<<defValue<<")";
   desc = descStream.str();
}//}}}
void ArgumentParser::usage(){//{{{
   map<string,Option>::iterator it;
   vector<string>::iterator itV;
   Option opt;
   message("\nUsage: %s ",(programName).c_str());
   sort(compulsory.begin(),compulsory.end());
   for(itV=compulsory.begin();itV!=compulsory.end();itV++){
      if(validOptions[*itV].shortName!="")
         message("-%s ",(validOptions[*itV].shortName).c_str());
      else
         message("--%s ",(validOptions[*itV].longName).c_str());
      if(validOptions[*itV].type!=OTBOOL)message("<%s> ",(*itV).c_str());
   }
   message(" [OPTIONS] %s\n",(argumentDesc).c_str());
   message("\n%s\n\nOptions:\n",(programDesc).c_str());
   message("  --help\n    Show this help information.\n\n");
   for(it=validOptions.begin();it!=validOptions.end();it++){
      opt=it->SS;
      message("  ");
      if(opt.shortName!=""){
         message("-%s",(opt.shortName).c_str());
         if(opt.type!=OTBOOL)message(" <%s>",(it->FF).c_str());
         if(opt.longName!="")message(" ,   ");
      }
      if(opt.longName!=""){
         message("--%s",(opt.longName).c_str());
         if(opt.type!=OTBOOL)message("=<%s>",(it->FF).c_str());
      }
      message("\n");
      if(opt.description!=""){
         message("    %s\n\n",(opt.description).c_str());
      }
   }
}//}}}
bool ArgumentParser::isSet(string name){//{{{
   if(! existsOption(name))return false;
   switch(validOptions[name].type){
      case OTSTRING:
         if(mapS.find(name)==mapS.end())return false;
         else return true;
      case OTLONG:
         if(mapL.find(name)==mapL.end())return false;
         else return true;
      case OTBOOL:
         if(mapB.find(name)==mapB.end())return false;
         else return mapB[name];
      case OTDOUBLE:
         if(mapD.find(name)==mapD.end())return false;
         else return true;
   }
   return false;
}//}}}
bool ArgumentParser::existsName(string name){//{{{
   if(names.find(name)==names.end())return false;
   return true;
}//}}}
bool ArgumentParser::existsOption(string name){//{{{
   if(validOptions.find(name)==validOptions.end())return false;
   return true;
}//}}}
