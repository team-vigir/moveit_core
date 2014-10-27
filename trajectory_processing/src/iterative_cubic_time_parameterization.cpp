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
        std::vector<double>        intervals;
        std::vector< std::vector<double> >      a_vec, b_vec, c_vec, d_vec, v_vec;
        intervals.resize(length-1, 0.0);
        a_vec.resize(length-2 , std::vector<double>(num_joints, 0.0));
        b_vec.resize(length-2 , std::vector<double>(num_joints, 0.0));
        c_vec.resize(length-2 , std::vector<double>(num_joints, 0.0));
        d_vec.resize(length-2 , std::vector<double>(num_joints, 0.0));
        v_vec.resize(length-2 , std::vector<double>(num_joints, 0.0));


        for (int i=0; i < (length -1); ++i)
        {
            int traj_ndx = start_ndx + i;
            intervals[i] = rob_trajectory.getWayPointDurationFromPrevious(traj_ndx+1);
        } //getVariablePosition

        for (int j=0; j< num_joints;++j)
        {
            int traj_ndx=0;
            logWarn("   Orig point %d  Joint %d Posn=%g Vel=%g Acc=%g  ",
                traj_ndx,j, rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(idx[j]),
                rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(idx[j]), rob_trajectory.getWayPointPtr(traj_ndx)->getVariableAcceleration(idx[j]));
        }

        for (int i=0; i < (length -2); ++i)
        {
            int traj_ndx = start_ndx + i;
            double time_factor_0 = (3.0/(intervals[i]*intervals[i+1]));
            double time_factor_1 = time_factor_0 * (intervals[i+1]*intervals[i+1]);
            time_factor_0 *= (intervals[i]*intervals[i]);

            for (int j=0; j < num_joints; ++j)
            {
                c_vec[i][j] = intervals[i];
                if (i < (length-3))
                {
                    a_vec[i+1][j] = intervals[i+2];
                }
                b_vec[i][j] = (2.0*(intervals[i] + intervals[i+1]));
                double q2 = rob_trajectory.getWayPointPtr(traj_ndx+2)->getVariablePosition(idx[j]);
                double q1 = rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariablePosition(idx[j]);
                double q0 = rob_trajectory.getWayPointPtr(traj_ndx  )->getVariablePosition(idx[j]);
                d_vec[i][j] = ( (q2 - q1)*time_factor_0 + (q1 - q0)*time_factor_1 );
                logWarn("   Orig point %d  Joint %d Posn=%g Vel=%g Acc=%g  (q2=%g,q1=%g,q0=%g,d_vec[i][j]=%g)",
                        i,j, rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariablePosition(idx[j]),
                        rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariableVelocity(idx[j]), rob_trajectory.getWayPointPtr(traj_ndx+1)->getVariableAcceleration(idx[j]),
                        q2,q1,q0,d_vec[i][j]);

            }
        }
        for (int j=0; j< num_joints; ++j)
        {
            d_vec[0][j]     -= rob_trajectory.getWayPointPtr(start_ndx)->getVariableVelocity(idx[j])*intervals[1];
            d_vec.back()[j] -= rob_trajectory.getWayPointPtr(length-1)->getVariableVelocity(idx[j])*intervals[length > 2 ? length-3 : 0];
            int traj_ndx=d_vec.size()-1;
            logWarn("   Orig point %d  Joint %d Posn=%g Vel=%g Acc=%g  ",
                    traj_ndx,j, rob_trajectory.getWayPointPtr(traj_ndx)->getVariablePosition(idx[j]),
                    rob_trajectory.getWayPointPtr(traj_ndx)->getVariableVelocity(idx[j]), rob_trajectory.getWayPointPtr(traj_ndx)->getVariableAcceleration(idx[j]));
        }

        // ----------------------------------------------------------
        // ------------------ Tri-diagonal solve -------------------
        int n = (int)d_vec.size();
        logWarn("  tridiagonal solve with n=%d", n);

        // forward elimination
        for (int i=1; i<n; i++)
        {
            for (int j=0; j < num_joints; ++j)
            {
               double m  = a_vec[i][j] / b_vec[i-1][j];
               b_vec[i][j] -= m*c_vec[i-1][j];
               d_vec[i][j] -= m*d_vec[i-1][j];
            }
        }

        // backward substitution
        for (int j=0; j < num_joints; ++j)
        {
            v_vec[n-1][j] = d_vec[n-1][j]/b_vec[n-1][j];
            logWarn("   v_vec[%d][%d] = %g  d=%g b=%g", n-1, j, v_vec[n-1][j], d_vec[n-1][j], b_vec[n-1][j]);

        }
        for (int i=n-2; i>=0; i--)
        {
            for (int j=0; j < num_joints; ++j)
            {
                v_vec[i][j] = (d_vec[i][j] - c_vec[i][j]*v_vec[i+1][j])/b_vec[i][j];
                logWarn("   v_vec[%d][%d] = %g  d=%g c=%g b=%g v[i+1][j]=%g", i, j, v_vec[i][j], d_vec[i][j], c_vec[i][j], b_vec[i][j], v_vec[i+1][j]);
            }
        }
        // ---- End tridiagonal solve for velocities

        // Assign smoothed velocities
        logWarn(" assign smoothed velocities with %d points", length);
        for (int i=0; i < (length -2); ++i)
        {
            int traj_ndx = start_ndx + i + 1;
            for (int j=0; j < num_joints; ++j)
            {
                logWarn(" assign smoothed velocities to traj_ndx= %d v_vec[%d][%d] = %g", traj_ndx,i,j,v_vec[i][j]);

                if (isnan(v_vec[i][j]))
                {
                    logError("    invalid v_vec[%d][%d] = %g", traj_ndx,i,j,v_vec[i][j]);
                    v_vec[i][j] = 0.00000123456;
                }
                rob_trajectory.getWayPointPtr(traj_ndx)->setVariableVelocity(idx[j],  v_vec[i][j]);

                // Calc coefficients for normalized time a t^3 + b t^2 + c t + d and solve for acceleration at final point

            }
        }

        start_ndx += length;
        if (start_ndx < num_points) --start_ndx;
        logWarn(" next start_ndx =  %d ", start_ndx);
    }
    logWarn("Completed velocity calculations! - now set accelerations");
    // Solve for the cubic spline coefficients
    for (int j=0; j < num_joints; ++j)
    {
        logWarn("   Point %d  Joint %d Posn=%g Vel=%g Acc=%g", 0,j, rob_trajectory.getWayPointPtr(0)->getVariablePosition(idx[j]), rob_trajectory.getWayPointPtr(0)->getVariableVelocity(idx[j]), rob_trajectory.getWayPointPtr(0)->getVariableAcceleration(idx[j]));
    }
    for (int i=1; i < num_points; ++i)
    {
        double dt = rob_trajectory.getWayPointDurationFromPrevious(i);
        double dt2 = dt*dt;

        for (int j=0; j < num_joints; ++j)
        {
            double d  = rob_trajectory.getWayPointPtr(i-1)->getVariablePosition(idx[j]);
            double c  = rob_trajectory.getWayPointPtr(i-1)->getVariableVelocity(idx[j])*dt;
            double b0 = rob_trajectory.getWayPointPtr(i)->getVariablePosition(idx[j]) - c - d;
            double b1 = rob_trajectory.getWayPointPtr(i)->getVariableVelocity(idx[j])*dt - c;
            double a  = -2.0*b0 + b1;
            double b  =  3.0*b0 - b1;

            if (i > 1)
            { // already specified acceleration from the previous segment
                if (fabs(2.0*b/dt2 - rob_trajectory.getWayPointPtr(i)->getVariableAcceleration(idx[j])) > 1e-6)
                {
                    logWarn("Inconsistent acceleration at point %d for joint %d  2b = %g   acc=%g ", i,j, 2*b/dt, rob_trajectory.getWayPointPtr(i)->getVariableAcceleration(idx[j]));
                }
            }
            else
            {
                rob_trajectory.getWayPointPtr(i-1)->setVariableAcceleration(idx[j], 2.0*b/dt2);
            }
            rob_trajectory.getWayPointPtr(i)->setVariableAcceleration(idx[j],  (6.0*a + 2.0*b)/dt2);
            logWarn("   Point %d  Joint %d Posn=%g Vel=%g Acc=%g  (a=%g,b=%g,c=%g,d=%g,b0=%g,b1=%g) dt=%g",
                    i,j, rob_trajectory.getWayPointPtr(i)->getVariablePosition(idx[j]), rob_trajectory.getWayPointPtr(i)->getVariableVelocity(idx[j]), rob_trajectory.getWayPointPtr(i)->getVariableAcceleration(idx[j]),
                    a,b,c,d,b0,b1,dt);

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
