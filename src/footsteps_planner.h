#pragma once
#include <ctime>
#include <chrono>
#include <mc_control/mc_controller.h>
#include "eigen-quadprog/eigen_quadprog_api.h"
#include "eigen-quadprog/QuadProg.h"
#include <lipm_walking/utils/polynomials.h>
#include <eigen3/Eigen/Dense>

namespace mc_plugin
{
namespace footsteps_planner
{

struct Admissible_Region {

    public :
        Admissible_Region();
        Admissible_Region(const Eigen::Vector3d & center,const Eigen::Vector3d & size){
            _center = center; _angle = _center.z(); _size = size; _size.z() = 0;
            _center.z() = 0;
            R.setZero();
            R(0,0) = cos(_angle) ; R(0,1) = -sin(_angle);
            R(1,0) = sin(_angle) ; R(1,1) =  cos(_angle);
            R(2,2) = 1;
            upper_left_corner  = _center + R * Eigen::Vector3d{-_size.x()/2, _size.y()/2,0};
            upper_right_corner = _center + R * Eigen::Vector3d{ _size.x()/2, _size.y()/2,0};
            lower_left_corner  = _center + R * Eigen::Vector3d{-_size.x()/2,-_size.y()/2,0};
            lower_right_corner = _center + R * Eigen::Vector3d{ _size.x()/2,-_size.y()/2,0};
            
            std::vector<Eigen::Vector3d> Polygone_Corners = {upper_left_corner,upper_right_corner,lower_right_corner,lower_left_corner};
            Polygone_Normals.resize(2,Polygone_Corners.size());
            Polygone_Edges_Center.resize(2,Polygone_Corners.size());
            Polygone_Vertices.resize(2,Polygone_Corners.size());
            Offset.resize(Polygone_Corners.size());
            for (int c = 0 ; c < Polygone_Corners.size() ; c++){
                Eigen::Vector3d point_1 = Polygone_Corners[c];
                Eigen::Vector3d point_2 = Polygone_Corners[ (c+1)%Polygone_Corners.size() ];
                Eigen::Vector3d vertice = (point_2 - point_1).normalized();
                Eigen::Vector3d normal = Eigen::Vector3d{0,0,1}.cross(vertice).normalized();
                Polygone_Normals(0,c) = normal.x();
                Polygone_Normals(1,c) = normal.y();
                Polygone_Vertices(0,c) = vertice.x();
                Polygone_Vertices(1,c) = vertice.y();


                Eigen::Matrix2d R_Vertices_0;
                R_Vertices_0(0,0) = Polygone_Normals(0,c); R_Vertices_0(1,0) = Polygone_Vertices(0,c);
                R_Vertices_0(1,0) = Polygone_Normals(1,c); R_Vertices_0(1,1) = Polygone_Vertices(1,c);

                Offset(c) = (R_Vertices_0.transpose() * Eigen::Vector2d{ point_1.x(),point_1.y()}).x();

 
            }

        }
        ~Admissible_Region() = default;

        std::vector<Eigen::Vector3d> Get_corners(){
            return {upper_left_corner,upper_right_corner,lower_right_corner,lower_left_corner};
        }

        Eigen::MatrixXd Polygone_Normals;
        Eigen::VectorXd Offset;

    private:
        
        Eigen::Vector3d _center;
        Eigen::Vector3d _size;
        double _angle;
        Eigen::Matrix3d R;
        Eigen::Vector3d upper_left_corner;
        Eigen::Vector3d upper_right_corner;
        Eigen::Vector3d lower_left_corner;
        Eigen::Vector3d lower_right_corner;
        Eigen::MatrixXd Polygone_Vertices;
        Eigen::MatrixXd Polygone_Edges_Center;

};

struct Steps_timings_output{

    public :
        Steps_timings_output(){};
        ~Steps_timings_output() = default;
        bool QPsuccess;
        Eigen::VectorXd Ts;
        double loss;

};

struct ref_traj_point{

    ref_traj_point(Eigen::Vector2d pose , double ori)
    {
        pose_ = pose;
        ori_ = ori;
    }
    ref_traj_point(sva::PTransformd pose)
    {
        pose_ = pose.translation().segment(0,2);
        ori_ = mc_rbdyn::rpyFromMat(pose.rotation()).z();
    }
    ~ref_traj_point() = default;

    Eigen::Vector2d pose()
    {
        return pose_;
    }
    Eigen::Vector3d vec3_pose()
    {
        return Eigen::Vector3d{pose_.x() , pose_.y() , 0.};
    }
    sva::PTransformd PT_pose() {
        Eigen::Vector3d center{pose_.x(),pose_.y(),0.};
        return sva::PTransformd(sva::RotZ(ori_) , center);
    }
    const double ori()
    {
        return ori_;
    }
    Eigen::Vector3d rpy_ori()
    {
        return Eigen::Vector3d{0. , 0. , ori_};
    }
    void ori(const double theta)
    {
        ori_ = theta;
    }

    Eigen::Vector2d pose_;
    double ori_;
};

struct Footstep : public mc_plugin::footsteps_planner::ref_traj_point{

    public:
        Footstep();
        Footstep(const sva::PTransformd & pose, double ts , const Eigen::Vector2d & step_size) : mc_plugin::footsteps_planner::ref_traj_point(pose)
        {
            ts_ = ts;
            Eigen::Vector3d center{pose_.x(),pose_.y(),ori_};
            Eigen::Vector3d dim{step_size_.x(),step_size_.y(),0.};
            step_rect_ = Admissible_Region(center,dim);
        }
        Footstep(const Eigen::Vector2d & pose,const double ori, double ts, const Eigen::Vector2d & step_size) : mc_plugin::footsteps_planner::ref_traj_point(pose,ori)
        {
            step_size_ = step_size;

            ts_ = ts;
            Eigen::Vector3d center{pose_.x(),pose_.y(),ori_};
            Eigen::Vector3d dim{step_size_.x(),step_size_.y(),0.};
            step_rect_ = Admissible_Region(center,dim);
        }
        ~Footstep() = default;

        Admissible_Region & rect() noexcept{
            return step_rect_;
        }
 
        const double ts() const noexcept   
        {
            return ts_;
        }

    private :

        Eigen::Vector2d pose_;
        Eigen::Vector2d step_size_;
        double ori_;
        Admissible_Region step_rect_;
        double ts_;

        
};

struct Footsteps_plan{

    public:
    
        Footsteps_plan();
        inline Footsteps_plan(const Footstep & support_foot , const Footstep & initial_swing_foot , const std::vector<Footstep> & footsteps)
        {
            footsteps_ = footsteps;
            support_foot_ = support_foot;
            initial_swing_foot_ = initial_swing_foot;       
        }
        ~Footsteps_plan() = default;

        void push_back(const Footstep & step)
        {
            footsteps_.push_back(step);
        }
        void clear()
        {
            footsteps_.clear();
        }
        std::vector<std::vector<Eigen::Vector3d>> get_steps_corners()
        {
         std::vector<std::vector<Eigen::Vector3d>> Output;
         for (int k  = 0 ; k < n_steps() ; k ++ )
         {
            Output.push_back(footsteps_[k].rect().Get_corners());
         }  
         return Output;
        }
        std::vector<Eigen::Vector3d> steps_pose()
        {
            std::vector<Eigen::Vector3d> output;
            for (int k = 0 ; k < footsteps_.size() ; k++ )
            {   
                
                output.push_back(footsteps_[k].vec3_pose());
            }
            return output;
        } 
        std::vector<sva::PTransformd> steps_PTpose()
        {
            std::vector<sva::PTransformd> output;
            for (int k = 0 ; k < footsteps_.size() ; k++ )
            {   
                
                output.push_back(footsteps_[k].PT_pose());
            }
            return output;
        } 
        std::vector<double> steps_timings()
        {
            std::vector<double> output;
            for (int k = 0 ; k < footsteps_.size() ; k++ )
            {   
                
                output.push_back(footsteps_[k].ts());
            }
            return output;
        }
        int n_steps() const noexcept
        {
            return footsteps_.size();
        }
        Footstep & support_foot()
        {
            return support_foot_;
        }
        Footstep & footstep(int indx)
        {
            return footsteps_[indx];
        }
        void support_foot(const Footstep & footstep)
        {
            support_foot_ = footstep;
        }

    
    private:

        std::vector<Footstep> footsteps_;
        Footstep support_foot_;
        Footstep initial_swing_foot_;


};

struct FootStepGen {
    
    public:
        FootStepGen(const mc_rtc::Configuration & config);
        FootStepGen();



        ~FootStepGen() = default;
        /**
         * Initialize the footsteps Generator
         * @tparam supportFootName 
         * @tparam P_f0 Support Foot
         * @tparam V Reference velocity inputs from t0 (can be empty)
         * @tparam Tstep desired ordered Steps Timing (can be empty)
         * @tparam Pf desired ordered Footsteps coordinate (can be empty) (angle in z coordinate)  
         */
        void Init(std::string supportFootName,
                    Footstep P_f0,
                    const std::vector<sva::MotionVecd> & V,
                    const std::vector<double> & Tstep,
                    std::vector<Footstep> & Pf);

        //Return The footsteps Theta values
        const Eigen::VectorXd & Theta_f() const noexcept{
            return Theta_f_;
        }    
    
        //Return the reference trajectory in the preview horizon
        std::vector<Eigen::Vector3d> Ref_Traj()
        {
            std::vector<Eigen::Vector3d> Output;
            for (int k = 0 ; k < P_traj_.size() ; k++)
            {
                Output.push_back(P_traj_[k].vec3_pose());
            }
            return Output;
        }
        //Compute The Footsteps and the Steps Timings
        Footsteps_plan compute_plan();

        Footsteps_plan & footsteps_plan()
        {
            return plan_;
        }
        
        const int & Get_Nsteps() const {
            return N_steps;
        }
        
    private:
  
        Eigen::VectorXd solveQP();

        /**
         * Compute N points trajectory between P_s_0 and P_s_1
         * @return Points Coordonate and angle of the trajectory
         */
        std::vector<ref_traj_point> GetRefTrajectory(ref_traj_point & P_s_0, ref_traj_point & P_s_1);

        std::vector< std::vector<double> > GetVelocityProfile(const Eigen::Vector3d & P_s_0, double V_Max, double V_Min, const std::vector<Eigen::Vector3d> & Traj);
        
        //Compute the Steps Timing dependings of the given parameter
        void GetStepsTimings();
        Steps_timings_output Get_constrained_Ts(const Eigen::VectorXd & Ts_candidate,
                                                const std::vector<Eigen::Vector2d> & StepsTimings_Upper_Lower_cstr);
        
        

        /**
         * return the position of the reference velocity integratin the velocity profile
         * @tparam k_end time index desired
         * @return Coordinate of the integrated ref velocity at time index k with orientation in z
         */
        ref_traj_point IntegrateVelProfile(int k_end); 
        int Get_ki(int k, int kfoot);
        
        std::string supportFoot = "RightFoot";

        HoubaPolynomial<Eigen::Vector2d> path;

        Footsteps_plan plan_;
        std::vector<Footstep> steps_inputs_;
        std::vector<double> t_steps_inputs_ ; //Input Step Timing
        std::vector<sva::MotionVecd> v_inputs_ ; //Velocity input

        int F_ ; //footsteps number   
        int N_steps = -1;

        std::vector<double> StepsTimings_; //Contains the time of each steps
        std::vector<int> StepsTimings_indx_; //Index timing of each steps
        std::vector<int> FootSteps_indx_; //Index of the input steps position for the right step timing
        Eigen::VectorXd Theta_f_ ; //Output Steps Angles
        Eigen::VectorXd m_Ts ;  //Steps Duration

        std::vector<ref_traj_point> P_traj_; //Position of reference trajectory for each timesteps 
                        
        //QP Problem
        bool QPsuccess = false;
        Eigen::MatrixXd Q_; //QP Hessian
        Eigen::VectorXd p_; //QP Grad

        Eigen::MatrixXd Aeq; //Equality Matrix
        Eigen::VectorXd beq; //Equality Vector

        Eigen::MatrixXd Aineq; //Inequality Matrix
        Eigen::VectorXd bineq; //Inequality Vector


        double Ts_min_ = 0.8 ; //Step Time lenght limits
        double Ts_max_ = 2; //Step Time lenght limits
        double l_ = 0.2; //Distance between foot
        double Tp_ = 3; // Preview horizon time
        double delta_ = 5e-2; //t_k - t_k-1
        double d_h_x = 0.2; //Next step tolerance zone
        double d_h_y = 0.05; //Next step tolerance zone
        double v_ = 0.1 ; //Cruise Parameters
        double max_theta = 3.14/6; //Max angle between two steps
        double P_; // Preview horizon time indexes
        double Ts_; //Cruise Parameters
        double robot_height_ = 150; //in cm

          
};

} //namespace foosteps_planner
} //namespace mc_plugin

