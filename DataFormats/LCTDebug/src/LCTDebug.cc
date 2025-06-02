#include "DataFormats/LCTDebug/interface/LCTDebug.h"


LCTDebugobject::LCTDebugobject(){
    bendingangle=-99;
    bendinganglenoalignmentcorrection = -99;
    KeyWG=0;
    bx=0;
    Bend=0;
    KeyStrip=0;
    slope=0;
}

void LCTDebugobject::Setbendingangle(int eighthStripDiff_){
    bendingangle = eighthStripDiff_;
}

void LCTDebugobject::Setbendinganglenoalignmentcorrection(int bendinganglenoalignmentcorrection_){
    bendinganglenoalignmentcorrection = bendinganglenoalignmentcorrection_;
}

void LCTDebugobject::Setidentifiers(int KeyWG_, int bx_, int Bend_,int KeyStrip_, int slope_){
    KeyWG = KeyWG_;
    bx = bx_;
    Bend = Bend_;
    KeyStrip = KeyStrip_;
    slope = slope_;
}

void LCTDebugobject::SetGEMClusterKeyStrip(int cl_es_){
    cl_es = cl_es_;
}

void LCTDebugobject::Setresidual(int residualwithalignment_){
    residualwithalignment = residualwithalignment_;
}

void LCTDebugobject::Setresidualnoalignmentcorrection(int residualwithoutalignment_){
    residualwithoutalignment = residualwithoutalignment_;
}

void LCTDebugobject::SetClusterRoll(int roll_){
    cluster_roll = roll_;
}

void LCTDebugobject::SetClusterBx(int bx_){
    cluster_bx = bx_;
}

int LCTDebugobject::Getbendingangle(){
    return bendingangle;
}

int LCTDebugobject::Getbendinganglenoalignmentcorrection(){
    return bendinganglenoalignmentcorrection;
}

std::vector<int> LCTDebugobject::Getidentifiers(){
    std::vector<int> output;
    output.push_back(KeyWG);
    output.push_back(bx);
    output.push_back(Bend);
    output.push_back(KeyStrip);
    output.push_back(slope);
    return output;
}

void LCTDebugobject::Setlayer2bool(bool isLayer2){
    layer2bool = isLayer2;
}

bool LCTDebugobject::Getlayer2bool(){
    return layer2bool;
}

int LCTDebugobject::GetGEMClusterKeyStrip(){
    return cl_es;
}

int LCTDebugobject::Getresidual(){
    return residualwithalignment;
}

int LCTDebugobject::Getresidualnoalignmentcorrection(){
    return residualwithoutalignment;
}

int LCTDebugobject::GetClusterRoll(){
    return cluster_roll;
}

int LCTDebugobject::GetClusterBx(){
    return cluster_bx;
}