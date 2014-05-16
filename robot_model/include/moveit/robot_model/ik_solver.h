/**************************************************************************
 * copyright 2014 SRI International.                                      *
 **************************************************************************/
/* Author: Dave Hershberger */

#ifndef MOVEIT_CORE_ROBOT_MODEL_IK_SOLVER_H
#define MOVEIT_CORE_ROBOT_MODEL_IK_SOLVER_H

namespace moveit
{
namespace core
{

/** @brief Superclass of inverse kinematics solver implementations.
 *
 * This is intended to supersede kinematics::KinematicsBase */
class IKSolver
{
public:
  typedef boost::function<bool(const RobotState&)> StateValidityCallbackFn;

  class Request
  {
  public:
    // seed_state_ is a pointer instead of anything fancier because I
    // want to make it easy for RobotState to send "this" in here,
    // without needing boost::shared_from_this and always creating
    // RobotStates in shared_ptrs.  That seemed pretty invasive.  Also
    // using a reference provides a bit more safety but ultimately is
    // not that different, and complicates initialization.

    /// State to take inspiration from.  Joints not being moved will
    /// be set to match the values in this.
    const RobotState* seed_state_;

    EigenSTL::vector_Affine3d target_poses_; // Goal poses for things, relative to reference_frame_.  Indexing matches getTipLinks().
    EigenSTL::vector_Affine3d targets_rel_tips_; // Poses of things relative to tip frames.  Indexing matches getTipLinks().
    std::string reference_frame_; // Name of coordinate frame targets are specified in.
    std::vector<std::string> lock_joints_; // Joints to avoid moving during this IK solution.

    // should this be a vector<double> instead?  Indexed same as group->getVariableNames?  or something?
    std::map<std::string, double> joint_difference_limits_; // how far from seed_state_ each joint value can be.

    StateValidityCallbackFn state_validity_callback_; // Function to accept or reject every solution.
    size_t num_attempts_; // Max number of restarts for iterative solvers. If 0 use default.
    double timeout_; // Max wall-clock time in seconds to spend on this call. If 0 use default.
    bool return_approximate_solution_; // ???
  };

  /** @brief Constructor.  Only default constructor is usable by
   * pluginlib::ClassLoader. */
  IKSolver() {}

  virtual ~IKSolver() {}

  /** @brief Implement this to perform group-specific initialization.
   * @param group The JointModelGroup this will do IK for.
   * @param config An XmlRpcValue of type XmlRpcValue::TypeStruct
   *        containing all the configuration parameters for this
   *        solver instance.  Typically the loader will read this from
   *        a rosparam corresponding to the group inside
   *        /robot_description_kinematics.
   * @param error_msg_out An optional place to put a descriptive error
   *        message when returning false (failure).
   * @return True on success, false on failure.
   *
   * Typically called during load of robot model.
   *
   * This should return false if the implementation will not work for
   * the current group and parameters, and should set *error_msg_out
   * to a description of the reason it won't work, if error_msg_out is
   * not NULL. */
  virtual bool initialize(const JointModelGroup* group,
                          const XmlRpc::XmlRpcValue& config,
                          std::string* error_msg_out = NULL) = 0;

  /** @brief Return the tip frames understood by this solver.
   *
   * When multiple tip frames are returned, the order is important.
   * The order of the tip frame names returned here is used to know
   * the meaning of Request::target_poses_ and
   * Request::targets_rel_tips_.
   *
   * These can be set by any means the implementation wants: read from
   * config data, computed based on the group, hard-coded, set at
   * run-time, whatever. */
  virtual std::vector<std::string> getTipFrames() const = 0;

  virtual void setDefaultTimeout(double timeout_seconds) { default_timeout_ = timeout_seconds; }
  virtual double getDefaultTimeout() const { return default_timeout_; }

  virtual void setDefaultAttempts(size_t num_attempts) { default_num_attempts_ = attempts; }
  virtual double getDefaultAttempts() const { return default_num_attempts_; }

  /** @brief Return the maximum number of simultaneous solutions this
   * solver can compute.
   *
   * Some IK solvers (like analytic solvers for 6DOF arms) can compute
   * more than one solution.  For instance, an "elbow-up" and an
   * "elbow-down" solution if the elbow joint can bend both ways.  The
   * value returned by this function should not depend on any
   * particular configuration (like being away from singularities).
   * It should always just return the maximum.
   *
   * If you don't know the exact number, it is safe to overestimate.
   * This is just available for users who want multiple solutions to
   * allocate the right number of RobotStates for their call to
   * solve().  If they allocate too many they will be safely
   * ignored.
   *
   * If you are writing an iterative solver the answer is probably
   * always 1, so that's what this default implementation returns. */
  virtual size_t getMaxPossibleSolutions() const { return 1; }

  /** @brief Try to solve IK for the given Request.
   * @param request contains the IK request.  It can be a pointer to a
   *        subclass of Request, in which case the implementation of
   *        IKSolver receiving it may or may not know about any extra
   *        information in it, depending on the implementation.
   * @param solutions_out provides the place where solve() will put
   *        valid solutions.  It will always fill it starting from
   *        index 0 and going up.  It will not attempt to change the
   *        number of elements in the vector.  The size of the vector
   *        passed indicates how many solutions (max) are desired by
   *        the caller.
   * @param error_msg_out optionally provides a place for solve() to
   *        put error messages.  Should only be used when returning 0.
   * @return the number of valid solutions that were actually found
   *         (and copied into solutions_out).  When returning 0, a
   *         descriptive error message should be copied into
   *         error_msg_out.
   *
   * Joint values for each solution are copied into each RobotState
   * pointed to from solutions_out.  If any of the RobotStatePtr s are
   * empty, the results are undefined.
   *
   * RobotState::update() is /em not called on the solution states.
   * It is up to the caller to do this if necessary. */
  virtual size_t solve(const RequestPtr& request,
                       const std::vector<RobotStatePtr>& solutions_out,
                       std::string* error_msg_out = NULL) const = 0;

  // Note: I used a vector of smart pointers instead of a vector of
  // RobotState instances in solve() so I wouldn't have to worry about
  // aligment of RobotState's Eigen member variables.
};

} // end namespace core
} // end namespace moveit

#endif // MOVEIT_CORE_ROBOT_MODEL_IK_SOLVER_H
