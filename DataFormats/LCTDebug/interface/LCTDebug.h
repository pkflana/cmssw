#ifndef LCTDEBUG_CLASSES
#define LCTDEBUG_CLASSES
#include <vector>

class LCTDebugobject{
    
    public:
        LCTDebugobject();
        void Setbendingangle(int eighthStripDiff_);
        void SetGEMClusterKeyStrip(int cl_es_);
        void Setresidual(int residualwithalignment);
        void Setresidualnoalignmentcorrection(int residualwithoutalignment_);
        void Setbendinganglenoalignmentcorrection(int bendinganglenoalignmentcorrection_);
        void Setidentifiers(int KeyWG_, int bx_, int Bend_,int KeyStrip_, int slope_);
        void SetClusterRoll(int roll_);
        void SetClusterBx(int bx_);
        int Getbendingangle();
        int Getbendinganglenoalignmentcorrection();
        void Setlayer2bool(bool isLayer2);
        bool Getlayer2bool();
        std::vector<int> Getidentifiers();
        int GetGEMClusterKeyStrip();
        int GetCLCTKeyStrip();
        int Getresidual();
        int Getresidualnoalignmentcorrection();
        int GetClusterRoll();
        int GetClusterBx();
    private:    
        int bendingangle;
        int bendinganglenoalignmentcorrection;
        int KeyWG;
        int bx;
        int Bend;
        int KeyStrip;
        int slope;
        bool layer2bool;
        int cl_es;
        int residualwithalignment;
        int residualwithoutalignment;
        int cluster_roll;
        int cluster_bx;
};











#endif