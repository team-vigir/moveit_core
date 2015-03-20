/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2009, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Authors: David Conner (conner@torcrobotics.com) based on work by Ken Anderson */

#include <moveit/trajectory_processing/iterative_cubic_time_parameterization.h>
#include <moveit_msgs/JointLimits.h>
#include <console_bridge/console.h>
#include <moveit/robot_state/conversions.h>

#include <Eigen/Dense>

namespace trajectory_processing
{

static const double DEFAULT_VEL_MAX = 1.0;
static const double DEFAULT_ACCEL_MAX = 1.0;
static const double ROUNDING_THRESHOLD = 0.01;

IterativeCubicTimeParameterization::IterativeCubicTimeParameterization(unsigned int max_iterations,
                                                                       double max_time_change_per_it)
  : IterativeParabolicTimeParameterization(max_iterations,max_time_change_per_it)
{}

IterativeCubicTimeParameterization::~IterativeCubicTimeParameterization()
{}

namespace
{

// Takes the time differences, and updates the timestamps, velocities and accelerations
// in the trajectory.
void updateTrajectory(robot_trajectory::RobotTrajectory& rob_trajectory,
                      const std::vector<double>& time_diff)
{
  // Error check
  if (time_diff.empty())
    return;

  double time_sum = 0.0;

  robot_state::RobotStatePtr prev_waypoint;
  robot_state::RobotStatePtr curr_waypoint;
  robot_state::RobotStatePtr next_waypoint;

  const robot_model::JointModelGroup *group = rob_trajectory.getGroup();
  const std::vector<std::string> &vars = group->getVariableNames();
  const std::vector<int> &idx = group->getVariableIndexList();

  int num_points = rob_trajectory.getWayPointCount();

  rob_trajectory.setWayPointDurationFromPrevious(0, time_sum);

  // Times
  for (int i = 1; i < num_points; ++i)
    // Update the time between the waypoints in the robot_trajectory.
    rob_trajectory.setWayPointDurationFromPrevious(i, time_diff[i-1]);

  // Return if there is only one point in the trajectory!
  if (num_points <= 1)
    return;

  // Accelerations
  for (int i = 0; i < num_points; ++i)
  {
    curr_waypoint = rob_trajectory.getWayPointPtr(i);

    if (i > 0)
      prev_waypoint = rob_trajectory.getWayPointPtr(i-1);

    if (i < num_points-1)
      next_waypoint = rob_trajectory.getWayPointPtr(i+1);

    for (std::size_t j = 0; j < vars.size(); ++j)
    {
      double q1;
      double q2;
      double q3;
      double dt1;
      double dt2;

      if (i == 0)
      {
        // First point
        q1 = next_waypoint->getVariablePosition(idx[j]);
        q2 = curr_waypoint->getVariablePosition(idx[j]);
        q3 = q1;

        dt1 = dt2 = time_diff[i];
      }
      else
        if (i < num_points-1)
        {
          // middle points
          q1 = prev_waypoint->getVariablePosition(idx[j]);
          q2 = curr_waypoint->getVariablePosition(idx[j]);
          q3 = next_waypoint->getVariablePosition(idx[j]);

          dt1 = time_diff[i-1];
          dt2 = time_diff[i];
        }
        else
        {
          // last point
          q1 = prev_waypoint->getVariablePosition(idx[j]);
          q2 = curr_waypoint->getVariablePosition(idx[j]);
          q3 = q1;

          dt1 = dt2 = time_diff[i-1];
        }

      double v1, v2, a;

      bool start_velocity = false;
      if (dt1 == 0.0 || dt2 == 0.0)
      {
        v1 = 0.0;
        v2 = 0.0;
        a = 0.0;
      }
      else
      {
        if (i == 0)
        {
          if (curr_waypoint->hasVelocities())
          {
            start_velocity = true;
            v1 = curr_waypoint->getVariableVelocity(idx[j]);
          }
        }
        v1 = start_velocity ? v1 : (q2-q1)/dt1;
        //v2 = (q3-q2)/dt2;
        v2 = start_velocity ? v1 : (q3-q2)/dt2; // Needed to ensure continuous velocity for first point
        a = 2.0*(v2-v1)/(dt1+dt2);
      }

      curr_waypoint->setVariableVelocity(idx[j], (v2+v1)/2.0);
      curr_waypoint->setVariableAcceleration(idx[j], a);
    }
  }
}
}

// This function recalculates the velocity and acceleration by assuming piecewise cubic splines between knot points.
// The function processes knot points in blocks of up to 45 knot points (stride).
// Two knots are added in each block (one after first point and one before last point) to balance the equations and unknowns.
// These extra knot points are not used in the final trajectory, which may introduce some perturbations  in interpolations.
void IterativeCubicTimeParameterization::smoothTrajectory(robot_trajectory::RobotTrajectory& rob_trajectory,
                                                          std::vector<double> & time_diff) const
{
    const int num_points = rob_trajectory.getWayPointCount();
    const robot_model::JointModelGroup *group = rob_trajectory.getGroup();
    const unsigned int num_joints = group->getVariableCount();
    const std::vector<int> &idx = group->getVariableIndexList();

    if (num_points <= 2)
    {
        logInform("     No cubic smoothing for trajectory with %u points!", num_points);
        return; // simple 2 point trajectory, no smoothing
    }

    int blocks = 1;
    int stride = num_points; // minimum length of smoothed points
    while (stride > 45)
    {
        ++blocks;
        stride = ceil(double(num_points)/double(blocks)) + 2; // allow room to back up during smoothing
    }

    logInform("   Do C3 smoothing of trajectory with %d points - use %d blocks with stride=%d ", num_points, blocks, stride );


    int start_ndx = 0;
    while (start_ndx < num_points)
    {
        int length    = num_points - start_ndx;

        if (length > stride)
        {
            length = stride;
        }
        int block_points = length + 2; // account for recalc of interior transition points
        logInform("     Process trajectory smoothing from %d to %d of %d total ", start_ndx, (start_ndx+length), num_points);

// p(tau) = a t^3 + b t^2 + c t + d
#define NP 4
#define a_ndx 0
#define b_ndx 1
#define c_ndx 2
#define d_ndx 3

        // Solve for parameters ===> M p = x ==> p = M^(-1)x
        int num_segments  = (block_points-1); // = length + 1
        int num_equations = NP*num_segments;
        Eigen::MatrixXd  M(num_equations, num_equations); M.setConstant(0.0);
        Eigen::VectorXd  x(num_equations); x.fill(0.0);


        // Get intervals (DT) for this block
        std::vector<double>        intervals;
        intervals.resize(num_segments,0.0);
        intervals[0] = 0.25*rob_trajectory.getWayPointDurationFromPrevious(start_ndx+1); // free transition point at start
        intervals[1] = 0.75*rob_trajectory.getWayPointDurationFromPrevious(start_ndx+1); // end of first segment
        ///double time_sum = intervals[0];
        for (int i=1; i < length-1; ++i)
        {
            int traj_ndx = start_ndx + i;
            intervals[i+1] = rob_trajectory.getWayPointDurationFromPrevious(traj_ndx+1);
            ///time_sum += intervals[i];
            ///logWarn("   interval[% 3d] = %g =?= %g  sum=%g",i,intervals[i],time_diff[i-1],time_sum);
        }
        intervals[length-1] *= 0.75;  // free transition point at end
        intervals[length] = 0.25*rob_trajectory.getWayPointDurationFromPrevious(start_ndx+length-1);
        ///time_sum += intervals[length-1];
        ///logWarn("   interval[% 3d] = %g =?= %g  sum=%g",length-1,intervals[length-1],0.75*time_diff[length-2],time_sum);
        ///time_sum += intervals[length];
        ///logWarn("   interval[% 3d] = %g =?= %g  sum=%g",length,  intervals[length]  ,0.25*time_diff[length-2],time_sum);


        // Fill in the coefficient matrix (M) that is the same for each joint
        M(0,d_ndx) = 1.0; //   q[0]
        M(1,c_ndx) = 1.0; //  dq[0]
        M(2,b_ndx) = 1.0; // ddq[0]
        int eqn_cnt = 3;

        // free transition point with continuous q,dq,ddq
        M(3,a_ndx) =     pow(intervals[0],3.0);  M(3,b_ndx) = pow(intervals[0],2.0); M(3,c_ndx) = intervals[0]; M(3,d_ndx) = 1.0; M(3,NP+d_ndx) = -1.0;  //continuous q at transition point
        M(4,a_ndx) = 3.0*pow(intervals[0],2.0);  M(4,b_ndx) = 2.0*intervals[0];      M(4,c_ndx) = 1.0;          M(4,d_ndx) = 0.0; M(4,NP+c_ndx) = -1.0;  //continuous dq at transition point
        M(5,a_ndx) = 6.0*    intervals[0];       M(5,b_ndx) = 2.0;                   M(5,c_ndx) = 0.0;          M(5,d_ndx) = 0.0; M(5,NP+b_ndx) = -2.0;  //continuous ddq at transition point
        eqn_cnt += 3;

        for (int i=1; i < num_segments-2; ++i)
        {
            // position at the end of current segment q(T)
            M(6+(i-1)*NP, i*NP + d_ndx) =         1.0;
            M(6+(i-1)*NP, i*NP + c_ndx) =     intervals[i];
            M(6+(i-1)*NP, i*NP + b_ndx) = pow(intervals[i],2.0);
            M(6+(i-1)*NP, i*NP + a_ndx) = pow(intervals[i],3.0);
            ++eqn_cnt;

            // position at start of next segment q(0)
            M(7+(i-1)*NP,(i+1)*NP + d_ndx) = 1.0;
            ++eqn_cnt;

            // continuous velocity at knot point dq(T)=dq(0)
            M(8+(i-1)*NP,    i*NP + c_ndx) =         1.0;
            M(8+(i-1)*NP,    i*NP + b_ndx) = 2.0 *   intervals[i];
            M(8+(i-1)*NP,    i*NP + a_ndx) = 3.0*pow(intervals[i],2.0);
            M(8+(i-1)*NP,(i+1)*NP + c_ndx) = -1.0;
            ++eqn_cnt;

            // continuous acceleration at knot point ddq(T)=ddq(0)
            M(9+(i-1)*NP,    i*NP + b_ndx) = 2.0 ;
            M(9+(i-1)*NP,    i*NP + a_ndx) = 6.0*intervals[i];
            M(9+(i-1)*NP,(i+1)*NP + b_ndx) = -2.0;
            ++eqn_cnt;

        }

        // penultimate knot point is free transition point with continuous q,dq,ddq
        int i = num_segments - 2;
        // free transibtion point with continuous q,dq,ddq
        M(6+(i-1)*NP,i*NP + a_ndx) =     pow(intervals[i],3.0);  M(6+(i-1)*NP,i*NP + b_ndx) = pow(intervals[i],2.0); M(6+(i-1)*NP,i*NP + c_ndx) = intervals[i]; M(6+(i-1)*NP,i*NP + d_ndx) = 1.0; M(6+(i-1)*NP,(i+1)*NP + d_ndx) = -1.0;  //continuous q at transition point
        M(7+(i-1)*NP,i*NP + a_ndx) = 3.0*pow(intervals[i],2.0);  M(7+(i-1)*NP,i*NP + b_ndx) = 2.0*intervals[i];      M(7+(i-1)*NP,i*NP + c_ndx) = 1.0;          M(7+(i-1)*NP,i*NP + d_ndx) = 0.0; M(7+(i-1)*NP,(i+1)*NP + c_ndx) = -1.0;  //continuous dq at transition point
        M(8+(i-1)*NP,i*NP + a_ndx) = 6.0*    intervals[i];       M(8+(i-1)*NP,i*NP + b_ndx) = 2.0;                   M(8+(i-1)*NP,i*NP + c_ndx) = 0.0;          M(8+(i-1)*NP,i*NP + d_ndx) = 0.0; M(8+(i-1)*NP,(i+1)*NP + b_ndx) = -2.0;  //continuous ddq at transition point
        eqn_cnt += 3;

        // Final knot point
        M( 9+(i-1)*NP,(i+1)*NP + a_ndx) =     pow(intervals[i],3.0);  M( 9+(i-1)*NP,(i+1)*NP + b_ndx) = pow(intervals[i],2.0); M( 9+(i-1)*NP,(i+1)*NP + c_ndx) = intervals[i]; M( 9+(i-1)*NP,(i+1)*NP + d_ndx) = 1.0; // q at transition point
        M(10+(i-1)*NP,(i+1)*NP + a_ndx) = 3.0*pow(intervals[i],2.0);  M(10+(i-1)*NP,(i+1)*NP + b_ndx) = 2.0*intervals[i];      M(10+(i-1)*NP,(i+1)*NP + c_ndx) = 1.0;          M(10+(i-1)*NP,(i+1)*NP + d_ndx) = 0.0; // dq at transition point
        M(11+(i-1)*NP,(i+1)*NP + a_ndx) = 6.0*    intervals[i];       M(11+(i-1)*NP,(i+1)*NP + b_ndx) = 2.0;                   M(11+(i-1)*NP,(i+1)*NP + c_ndx) = 0.0;          M(11+(i-1)*NP,(i+1)*NP + d_ndx) = 0.0; // ddq at transition point
        eqn_cnt += 3;

        if (eqn_cnt != num_equations)
        {
            logError("      Invalid number of equations %d =/= %d",eqn_cnt,num_equations);
        }

        if ((11+(i-1)*NP) != (num_equations-1))
        {
            logError("      Invalid indexing last row = %d =/= number of equations-1= %d ",(11+(i-1)*NP),(num_equations-1));

        }
        //std::cout << "M=" << M << std::endl;

        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Qr = M.colPivHouseholderQr();
        if (!Qr.isInvertible())
        {
            //std::cout << "M=" << M << std::endl;
            logError("       Matrix Qr  rank=%d invertible=%d - smoothing failed at start_ndx=%d of %d total points!",
                    Qr.rank(), Qr.isInvertible(),start_ndx,num_points);
            logError("       num_equations=%d eqn_cnt = %d number segments=%d ",num_equations, eqn_cnt, num_segments);
            std::cerr <<" intervals=[";
            for (int i=0; i < intervals.size(); ++i)
            {
                std::cerr << intervals[i];
            }
            std::cerr << "]" << std::endl;
            return;
        }

        for (int j=0; j< num_joints;++j)
        {
            int jnt = idx[j];

            ///logWarn("       Process joint %d (%d)...",j,jnt);
            // Now populate the "x" vector of Mp=x for each joint, and solve for coefficients of that joint
            x.fill(0.0); // clear for each joint
            x(0) = rob_trajectory.getWayPointPtr(start_ndx)->getVariablePosition(jnt);
            x(1) = rob_trajectory.getWayPointPtr(start_ndx)->getVariableVelocity(jnt);
            x(2) = rob_trajectory.getWayPointPtr(start_ndx)->getVariableAcceleration(jnt);

            // x(3)...x(5) are continuous (delta=0)
            eqn_cnt = 6;

            int i, traj_ndx;
            for (i=1; i < num_segments-2; ++i)
            {
                traj_ndx = start_ndx + i;

                // Defined positions on interior points
                x(6+(i-1)*NP) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(jnt);
                x(7+(i-1)*NP) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(jnt);

                // x(8+), x(9+) are continuous dq and ddq  ==> delta = 0.0

                eqn_cnt += 4;
            }

            i = num_segments - 2; // = length -1

            // penultimate point (x(6)..x(8)) has 3 continuity eqns delta = 0.0
            eqn_cnt += 3;

            // Final point
            traj_ndx = start_ndx + i;
            //logWarn("       set final block point data with i=%d traj_ndx=%d",i,traj_ndx);
            x( 9+(i-1)*NP) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(jnt);

            // Velocity at final segment
            if (traj_ndx == (num_points-1))
            {
                //logWarn("           set terminal point data with i=%d traj_ndx=%d",i,traj_ndx);
                x(10+(i-1)*NP) = 0.0;  // assume we always come to rest at final point
                x(11+(i-1)*NP) = 0.0;  // assume we always come to rest at final point
            }
            else
            {
                x(10+(i-1)*NP) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(jnt);
                x(11+(i-1)*NP) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariableAcceleration(jnt);
            }
            eqn_cnt += 3;

            if (eqn_cnt != num_equations)
            {
                logError(" Invalid number of equations %d =/= %d for joint %d (last row=%d)",eqn_cnt,num_equations,j,(11+(i-1)*NP));
            }

            // Calc coefficients for segment parameters using a t^3 + b t^2 + c t + d
            Eigen::VectorXd parameters = Qr.solve(x);

            // Update velocity and acceleration date for the interior knot points
            for (i=1; i< length-1; ++i)
            {
                traj_ndx = start_ndx + i; // next point in trajectory
                int seg = i + 1;

                ///double a = parameters(seg * NP + a_ndx); //a t^3 + b t^2 + c t + d
                double b = parameters(seg * NP + b_ndx);
                double c = parameters(seg * NP + c_ndx);
                ///double d = parameters(seg * NP + d_ndx);
                ///if (0==j)
                ///{
                ///    logWarn("     i=%d traj_ndx = %d joint=%d", i,traj_ndx, j);
                ///    logWarn("           a=%g b=%g c=%g d=%g =?= %g = q dt=%f",a,b,c,d, rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(jnt),intervals[i]);
                ///}
                // Set the interior knot values based on dq(0) and ddq(0) values
                rob_trajectory.getWayPointPtr(traj_ndx)->setVariableVelocity(    jnt,  c   );
                rob_trajectory.getWayPointPtr(traj_ndx)->setVariableAcceleration(jnt, 2.0*b);
            }            
            // Terminal data stays the same
        }
        start_ndx += length;
        ///logWarn("     new start = %d of %d total points (length=%d)",start_ndx, num_points, length);
        if (start_ndx != num_points)
        {
            // Re calculate starting from interior knot point of prior block
            start_ndx -= 2;
            logInform("  Calculate the next block starting at index = %d of %d total points (length=%d)",start_ndx, num_points, length);
        }
    }

}

bool IterativeCubicTimeParameterization::computeTimeStamps(robot_trajectory::RobotTrajectory& trajectory,
                                                           const double max_velocity_scaling_factor) const
{
  if (trajectory.empty())
    return true;

  const robot_model::JointModelGroup *group = trajectory.getGroup();
  if (!group)
  {
    logError("It looks like the planner did not set the group the plan was computed for");
    return false;
  }

  const int num_points = trajectory.getWayPointCount();
  const unsigned int num_joints = group->getVariableCount();

  // this lib does not actually work properly when angles wrap around, so we need to unwind the path first
  trajectory.unwind();

  std::vector<double> time_diff(num_points-1, 0.00005);   // the time difference between adjacent points - assign minimal value

  // Use iterative parabolic time parameterization to calculate the time differences given constraints and time scaling
  applyVelocityConstraints(trajectory, time_diff, max_velocity_scaling_factor);
  applyAccelerationConstraints(trajectory, time_diff);
  updateTrajectory(trajectory, time_diff);

  // Recalculate the velocity and accelerations based on
  // piecewise cubic splines with continuous q,dq,ddq at knot points
  smoothTrajectory(trajectory, time_diff);

  return true;
}

}
