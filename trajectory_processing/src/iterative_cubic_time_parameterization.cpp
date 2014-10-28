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

/* Author: Ken Anderson */

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
  : max_iterations_(max_iterations),
    max_time_change_per_it_(max_time_change_per_it)
{}

IterativeCubicTimeParameterization::~IterativeCubicTimeParameterization()
{}

namespace
{
void printPoint(const trajectory_msgs::JointTrajectoryPoint& point, std::size_t i) 
{
  logDebug(  " time   [%i]= %f",i,point.time_from_start.toSec());
  if(point.positions.size() >= 7 )
  {
    logDebug(" pos_   [%i]= %f %f %f %f %f %f %f",i,
             point.positions[0],point.positions[1],point.positions[2],point.positions[3],point.positions[4],point.positions[5],point.positions[6]);
  }
  if(point.velocities.size() >= 7 )
  {
    logDebug("  vel_  [%i]= %f %f %f %f %f %f %f",i,
             point.velocities[0],point.velocities[1],point.velocities[2],point.velocities[3],point.velocities[4],point.velocities[5],point.velocities[6]);
  }
  if(point.accelerations.size() >= 7 )
  {
    logDebug("   acc_ [%i]= %f %f %f %f %f %f %f",i,
             point.accelerations[0],point.accelerations[1],point.accelerations[2],point.accelerations[3],point.accelerations[4],point.accelerations[5],point.accelerations[6]);
  }
}

void printStats(const trajectory_msgs::JointTrajectory& trajectory, const std::vector<moveit_msgs::JointLimits>& limits)
{
  logDebug("jointNames= %s %s %s %s %s %s %s",
           limits[0].joint_name.c_str(),limits[1].joint_name.c_str(),limits[2].joint_name.c_str(),
           limits[3].joint_name.c_str(),limits[4].joint_name.c_str(),limits[5].joint_name.c_str(),
           limits[6].joint_name.c_str());
  logDebug("maxVelocities= %f %f %f %f %f %f %f",
           limits[0].max_velocity,limits[1].max_velocity,limits[2].max_velocity,
           limits[3].max_velocity,limits[4].max_velocity,limits[5].max_velocity,
           limits[6].max_velocity);
  logDebug("maxAccelerations= %f %f %f %f %f %f %f",
           limits[0].max_acceleration,limits[1].max_acceleration,limits[2].max_acceleration,
           limits[3].max_acceleration,limits[4].max_acceleration,limits[5].max_acceleration,
           limits[6].max_acceleration);
  // for every point in time:
  for (std::size_t i = 0; i<trajectory.points.size(); ++i)
    printPoint(trajectory.points[i], i);
}
}

// Applies velocity
void IterativeCubicTimeParameterization::applyVelocityConstraints(robot_trajectory::RobotTrajectory& rob_trajectory,
                                                                      std::vector<double> &time_diff) const
{
  const robot_model::JointModelGroup *group = rob_trajectory.getGroup();
  const std::vector<std::string> &vars = group->getVariableNames();
  const std::vector<int> &idx = group->getVariableIndexList();
  const robot_model::RobotModel &rmodel = group->getParentModel();
  const int num_points = rob_trajectory.getWayPointCount();

  for (int i = 0 ; i < num_points-1 ; ++i)
  {
    const robot_state::RobotStatePtr &curr_waypoint = rob_trajectory.getWayPointPtr(i);
    const robot_state::RobotStatePtr &next_waypoint = rob_trajectory.getWayPointPtr(i+1);
    
    for (std::size_t j = 0 ; j < vars.size() ; ++j)
    {
      double v_max = 1.0;
      const robot_model::VariableBounds &b = rmodel.getVariableBounds(vars[j]);
      if (b.velocity_bounded_)
        v_max = std::min(fabs(b.max_velocity_), fabs(b.min_velocity_));
      const double dq1 = curr_waypoint->getVariablePosition(idx[j]);
      const double dq2 = next_waypoint->getVariablePosition(idx[j]);
      const double t_min = std::abs(dq2-dq1) / v_max;
      if (t_min > time_diff[i])
        time_diff[i] = t_min;
    }
  }
}

// Iteratively expand dt1 interval by a constant factor until within acceleration constraint
// In the future we may want to solve to quadratic equation to get the exact timing interval.
// To do this, use the CubicTrajectory::quadSolve() function in cubic_trajectory.h
double IterativeCubicTimeParameterization::findT1(const double dq1,
                                                      const double dq2,
                                                      double dt1,
                                                      const double dt2,
                                                      const double a_max) const
{
  const double mult_factor = 1.01;
  double v1 = (dq1)/dt1;
  double v2 = (dq2)/dt2;
  double a = 2.0*(v2-v1)/(dt1+dt2);

  while( std::abs( a ) > a_max )
  {
    v1 = (dq1)/dt1;
    v2 = (dq2)/dt2;
    a = 2.0*(v2-v1)/(dt1+dt2);
    dt1 *= mult_factor;
  }

  return dt1;
}

double IterativeCubicTimeParameterization::findT2(const double dq1,
                                                      const double dq2,
                                                      const double dt1,
                                                      double dt2,
                                                      const double a_max) const
{
  const double mult_factor = 1.01;
  double v1 = (dq1)/dt1;
  double v2 = (dq2)/dt2;
  double a = 2.0*(v2-v1)/(dt1+dt2);

  while( std::abs( a ) > a_max )
  {
    v1 = (dq1)/dt1;
    v2 = (dq2)/dt2;
    a = 2.0*(v2-v1)/(dt1+dt2);
    dt2 *= mult_factor;
  }

  return dt2;
}

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


// Applies Acceleration constraints
void IterativeCubicTimeParameterization::applyAccelerationConstraints(robot_trajectory::RobotTrajectory& rob_trajectory,
                                                                          std::vector<double> & time_diff) const
{
  robot_state::RobotStatePtr prev_waypoint;
  robot_state::RobotStatePtr curr_waypoint;
  robot_state::RobotStatePtr next_waypoint;

  const robot_model::JointModelGroup *group = rob_trajectory.getGroup();
  const std::vector<std::string> &vars = group->getVariableNames();
  const std::vector<int> &idx = group->getVariableIndexList();
  const robot_model::RobotModel &rmodel = group->getParentModel();
  
  const int num_points = rob_trajectory.getWayPointCount();
  const unsigned int num_joints = group->getVariableCount();
  int num_updates = 0;
  int iteration = 0;
  bool backwards = false;
  double q1;
  double q2;
  double q3;
  double dt1;
  double dt2;
  double v1;
  double v2;
  double a;

  do
  {
    num_updates = 0;
    iteration++;

    // In this case we iterate through the joints on the outer loop.
    // This is so that any time interval increases have a chance to get propogated through the trajectory
    for (unsigned int j = 0; j < num_joints ; ++j)
    {
      // Loop forwards, then backwards
      for (int count = 0; count < 2; ++count)
      {
        for (int i = 0 ; i < num_points-1; ++i)
        {
          int index = backwards ? (num_points-1)-i : i;
          
          curr_waypoint = rob_trajectory.getWayPointPtr(index);

          if (index > 0)
            prev_waypoint = rob_trajectory.getWayPointPtr(index-1);
          
          if (index < num_points-1)
            next_waypoint = rob_trajectory.getWayPointPtr(index+1);

          // Get acceleration limits
          double a_max = 1.0;
          const robot_model::VariableBounds &b = rmodel.getVariableBounds(vars[j]);
          if (b.acceleration_bounded_)
            a_max = std::min(fabs(b.max_acceleration_), fabs(b.min_acceleration_));
          
          if (index == 0)
          { 
            // First point
            q1 = next_waypoint->getVariablePosition(idx[j]);
            q2 = curr_waypoint->getVariablePosition(idx[j]);
            q3 = next_waypoint->getVariablePosition(idx[j]);

            dt1 = dt2 = time_diff[index];
            assert(!backwards);
          }
          else
            if (index < num_points-1)
            {
              // middle points
              q1 = prev_waypoint->getVariablePosition(idx[j]);
              q2 = curr_waypoint->getVariablePosition(idx[j]);
              q3 = next_waypoint->getVariablePosition(idx[j]);
              
              dt1 = time_diff[index-1];
              dt2 = time_diff[index];
            }
            else
            { 
              // last point - careful, there are only numpoints-1 time intervals
              q1 = prev_waypoint->getVariablePosition(idx[j]);
              q2 = curr_waypoint->getVariablePosition(idx[j]);
              q3 = prev_waypoint->getVariablePosition(idx[j]);
              
              dt1 = dt2 = time_diff[index-1];
              assert(backwards);
            }
          
          if (dt1 == 0.0 || dt2 == 0.0)
          {
            v1 = 0.0;
            v2 = 0.0;
            a = 0.0;
          } 
          else
          {
            bool start_velocity = false;
            if (index == 0)
            {
              if (curr_waypoint->hasVelocities())
              {
                start_velocity = true;
                v1 = curr_waypoint->getVariableVelocity(idx[j]);
              }
            }
            v1 = start_velocity ? v1 : (q2-q1)/dt1;
            v2 = (q3-q2)/dt2;
            a = 2.0*(v2-v1)/(dt1+dt2);
          }

          if (fabs(a) > a_max + ROUNDING_THRESHOLD)
          {
            if (!backwards)
            {
              dt2 = std::min(dt2 + max_time_change_per_it_, findT2(q2-q1, q3-q2, dt1, dt2, a_max));
              time_diff[index] = dt2;
            }
            else
            {
              dt1 = std::min(dt1 + max_time_change_per_it_, findT1(q2-q1, q3-q2, dt1, dt2, a_max));
              time_diff[index-1] = dt1;
            }
            num_updates++;
            
            if (dt1 == 0.0 || dt2 == 0.0)
            {
              v1 = 0.0;
              v2 = 0.0;
              a = 0.0;
            } 
            else
            {
              v1 = (q2-q1)/dt1;
              v2 = (q3-q2)/dt2;
              a = 2*(v2-v1)/(dt1+dt2);
            }
          }
        }
        backwards = !backwards;
      }
    }
    //logDebug("applyAcceleration: num_updates=%i", num_updates);
  } while (num_updates > 0 && iteration < static_cast<int>(max_iterations_));
}


 // Based on code from http://docs.ros.org/fuerte/api/spline_smoother/html/clamped__cubic__spline__smoother_8h_source.html
void IterativeCubicTimeParameterization::smoothTrajectory(robot_trajectory::RobotTrajectory& rob_trajectory,
                                                          std::vector<double> & time_diff) const
{
    const int num_points = rob_trajectory.getWayPointCount();
    const robot_model::JointModelGroup *group = rob_trajectory.getGroup();
    const unsigned int num_joints = group->getVariableCount();
    const std::vector<int> &idx = group->getVariableIndexList();

    if (num_points <= 2)
    {
        return; // simple 2 point trajectory, no smoothing
    }

    logWarn("Smoothing trajectory with %d points",num_points);
    int start_ndx = 0;

    while (start_ndx < num_points)
    {
        int length    = num_points - start_ndx;

        if (length > 20)
        {
            length = 20; // limit to 20 points per http://docs.ros.org/fuerte/api/spline_smoother/html/clamped__cubic__spline__smoother_8h_source.html
            logWarn("Process trajectory smoothing from %d to %d in blocks of 20 points - will mess up 20th interior point", start_ndx, (start_ndx+length));
        }

// p(tau) = d tau^3 + c tau^2 + b tau + a,   where tau = dt/DT
#define NUM_PARAMETERS 4
#define a_ndx 0
#define b_ndx 1
#define c_ndx 2
#define d_ndx 3

        int num_segments  = (length-1);
        int num_equations = NUM_PARAMETERS*num_segments;
        Eigen::MatrixXd  A(num_equations, num_equations); A.setConstant(0.0);
        Eigen::VectorXd  b(num_equations); b.fill(0.0);


        // Get intervals (DT) for this block
        std::vector<double>        intervals;
        std::vector<double>        inv_dt;
        intervals.resize(num_segments,0.0);
        inv_dt.resize(num_segments,1.0);
        for (int i=0; i < num_segments; ++i)
        {
            int traj_ndx = start_ndx + i;
            intervals[i] = rob_trajectory.getWayPointDurationFromPrevious(traj_ndx+1);
            if (intervals[i] > 0.0001)
            {
                inv_dt[i] = 1.0/intervals[i];
            }
            else
            {
                logError("Invalid segment interval = %f",intervals[i]);
            }
        }

        // Fill in the coefficient matrix (A) that is the same for each joint
        int eqn_cnt = 0;
        for (int i=0; i < num_segments; ++i)
        {
            // Fill the the start position(tau=0)
            A(i,             i*NUM_PARAMETERS + a_ndx) = 1; // p(0)
            ++eqn_cnt;

            // Fill the end position (tau=1)
            A(i+num_segments,i*NUM_PARAMETERS + a_ndx) = 1; // p(1)
            A(i+num_segments,i*NUM_PARAMETERS + b_ndx) = 1;
            A(i+num_segments,i*NUM_PARAMETERS + c_ndx) = 1;
            A(i+num_segments,i*NUM_PARAMETERS + d_ndx) = 1;
            ++eqn_cnt;
        }
        // Fill the initial velocity
        int i = 0;
        A(i+2*num_segments,  i * NUM_PARAMETERS + a_ndx) = 0.0; // v(1))
        A(i+2*num_segments,  i * NUM_PARAMETERS + b_ndx) = inv_dt[i];
        A(i+2*num_segments,  i * NUM_PARAMETERS + c_ndx) = 0.0;
        A(i+2*num_segments,  i * NUM_PARAMETERS + d_ndx) = 0.0;
        ++eqn_cnt;

        // Velocity at final segment
        A(i+3*num_segments,  i * NUM_PARAMETERS + a_ndx) = 0.0; // v(1))
        A(i+3*num_segments,  i * NUM_PARAMETERS + b_ndx) = inv_dt.back();
        A(i+3*num_segments,  i * NUM_PARAMETERS + c_ndx) = 2.0*inv_dt.back();
        A(i+3*num_segments,  i * NUM_PARAMETERS + d_ndx) = 3.0*inv_dt.back();
        ++eqn_cnt;

        for (int i=1; i < (num_segments-1); ++i)
        {
            // Fill the equivalent velocities at interior point
            A(i+2*num_segments,  i * NUM_PARAMETERS + a_ndx) = 0.0; // v(1))
            A(i+2*num_segments,  i * NUM_PARAMETERS + b_ndx) = inv_dt[i-1];
            A(i+2*num_segments,  i * NUM_PARAMETERS + c_ndx) = 2.0*inv_dt[i-1];
            A(i+2*num_segments,  i * NUM_PARAMETERS + d_ndx) = 3.0*inv_dt[i-1];
            A(i+2*num_segments,(i+1)*NUM_PARAMETERS + a_ndx) =  0.0; // - v(0)
            A(i+2*num_segments,(i+1)*NUM_PARAMETERS + b_ndx) = -inv_dt[i];
            A(i+2*num_segments,(i+1)*NUM_PARAMETERS + c_ndx) =  0.0;
            A(i+2*num_segments,(i+1)*NUM_PARAMETERS + d_ndx) =  0.0;
            ++eqn_cnt;

            // Fill the equivalent accelerations at interior point
            A(i+3*num_segments,  i * NUM_PARAMETERS + a_ndx) = 0.0; // p(0)
            A(i+3*num_segments,  i * NUM_PARAMETERS + b_ndx) = 0.0;
            A(i+3*num_segments,  i * NUM_PARAMETERS + c_ndx) = 2.0*inv_dt[i-1];
            A(i+3*num_segments,  i * NUM_PARAMETERS + d_ndx) = 6.0*inv_dt[i-1];
            A(i+3*num_segments,(i+1)*NUM_PARAMETERS + a_ndx) =  0.0; // p(0)
            A(i+3*num_segments,(i+1)*NUM_PARAMETERS + b_ndx) =  0.0;
            A(i+3*num_segments,(i+1)*NUM_PARAMETERS + c_ndx) = -2.0*inv_dt[i];
            A(i+3*num_segments,(i+1)*NUM_PARAMETERS + d_ndx) =  0.0;
            ++eqn_cnt;
        }

        if (eqn_cnt != num_equations)
        {
            logError(" Invalid number of equations %d =/= %d",eqn_cnt,num_equations);
        }

        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Qr = A.colPivHouseholderQr();
        logWarn("  Matrix Qr  rank=%d invertible=%d",Qr.rank(), Qr.isInvertible());


        for (int j=0; j< num_joints;++j)
        {
            int jnt = idx[j];

            // Now populate the "b" vector of Ax=b for each joint, and solve for coefficients of that joint
            int i, traj_ndx, eqn_cnt = 0;
            for (i=0; i < num_segments; ++i)
            {
                traj_ndx = start_ndx + i;

                // Fill the the start position(tau=0)
                b(i) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(jnt);
                ++eqn_cnt;

                // Fill the end position (tau=1)
                b(i+num_segments) = rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariablePosition(jnt);
                ++eqn_cnt;
            }
            // Fill the initial velocity
            i = 0;
            traj_ndx = start_ndx + i;
            b(i+2*num_segments) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(jnt);
            ++eqn_cnt;

            // Velocity at final segment
            traj_ndx = start_ndx + num_segments;
            b(i+3*num_segments) = rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(jnt);
            ++eqn_cnt;

            // Differences between interior velocities and accelerations are zero (set with construction)
            eqn_cnt += (num_segments-1)*2;

            if (eqn_cnt != num_equations)
            {
                logError(" Invalid number of equations %d =/= %d",eqn_cnt,num_equations);
            }

            // Calc coefficients for normalized time d t^3 + c t^2 + b t + a
            Eigen::VectorXd parameters = Qr.solve(b);

            for (i=0; i< num_segments;++i)
            {

                traj_ndx = start_ndx + i;

                double a = parameters(i * NUM_PARAMETERS + a_ndx);
                double b = parameters(i * NUM_PARAMETERS + b_ndx);
                double c = parameters(i * NUM_PARAMETERS + c_ndx);
                double d = parameters(i * NUM_PARAMETERS + d_ndx);

                // Caclulate end point values
                double p0 = a;
                double p1 = d+c+b+a;
                double v0 = b*inv_dt[i];
                double v1 = (3.0*d+2.0*c+b)*inv_dt[i];
                double a0 = c;
                double a1 = (6.0*d + 2.0*c)*inv_dt[i]*inv_dt[i];

                // Sanity check positions
                double dp0 = p0 - rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(jnt);
                double dp1 = p1 - rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariableVelocity(jnt);
                if (fabs(dp0) > 1e-4 || fabs(dp1) > 1e-4)
                {
                    logError(" joint=%d segment %d joint position error dp0=%f dp1=%f", j, i, dp0, dp1);
                }

                // Calc velocity
                if (i)
                {
                    // check the prior starting values
                    double dv0 = v0 - rob_trajectory.getWayPointPtr(traj_ndx-1)->getVariableVelocity(jnt);
                    double da0 = a0 - rob_trajectory.getWayPointPtr(traj_ndx-1)->getVariableAcceleration(jnt);
                    if (fabs(dv0) > 1e-4*(v0 > 0 ? v0 : 1e-2) || fabs(da0) > 1e-4*(a0 > 0 ? a0 : 1e-2))
                    {
                        logError(" joint=%d segment %d joint velocity error dv0=%f acceleration error da0=%f", j, i, dv0, da0);
                    }
                }
                else
                {
                    // Set the initial values
                    double dv0 = v0 - rob_trajectory.getWayPointPtr(traj_ndx-1)->getVariableVelocity(jnt);
                    if (fabs(dv0) > 1e-4*(v0 > 0 ? v0 : 1e-2) )
                    {
                        logError(" joint=%d segment %d joint velocity error dv0=%f ", j, i, dv0);
                    }
                    rob_trajectory.getWayPointPtr(traj_ndx+1)->setVariableAcceleration(jnt,a0);
                }

                // Set the terminal values
                rob_trajectory.getWayPointPtr(traj_ndx+1)->setVariableVelocity(    jnt, v1);
                rob_trajectory.getWayPointPtr(traj_ndx+1)->setVariableAcceleration(jnt, a1);

            }
        }
    }
    logWarn("Smoothed trajectory with %d points and %d joints",num_points, num_joints);

}


bool IterativeCubicTimeParameterization::computeTimeStamps(robot_trajectory::RobotTrajectory& trajectory) const
{
  if (trajectory.empty())
    return true;

  const robot_model::JointModelGroup *group = trajectory.getGroup();
  if (!group)
  {
    logError("It looks like the planner did not set the group the plan was computed for");
    return false;
  }

  // this lib does not actually work properly when angles wrap around, so we need to unwind the path first
  trajectory.unwind();

  const int num_points = trajectory.getWayPointCount();
  std::vector<double> time_diff(num_points-1, 0.0);       // the time difference between adjacent points
  
  applyVelocityConstraints(trajectory, time_diff);
  applyAccelerationConstraints(trajectory, time_diff);
  
  updateTrajectory(trajectory, time_diff);
  smoothTrajectory(trajectory, time_diff);
  return true;
}

}
