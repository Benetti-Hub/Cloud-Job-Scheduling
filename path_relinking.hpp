#pragma once
#include "random_greedy.hpp"

//Here we define the type of search of the random greedy using the path relinking approach.
typedef std::list<std::pair<Job, int>> pool;

class Path_relinking : public Random_greedy{
  private:

  //List of parameters related to how PR is performed (Will determine the speed and the efficiency of the relinking process)
  unsigned K = 10;                    //pool size of the elite solutions;
  unsigned MAX_DEPTH = nodes.size();  //maximum number of steps in the relinking
  bool     FUTURE_SIGHT = true;       //explore one step further each time

  // parameters for randomization
  double alpha = 0.3;                  // parameter for the random setup selection
  double pi = 0.05;                    // parameter for random swaps in the queue
  
  pool moves = {};                     //Pool of moves

  //Map of the already explored move (hash function was made and is in job.hpp)
  std::unordered_map<std::pair<Job, int>, double> explored_moves = {};

  job_schedule_t best_relinking_sch;        //Best relinking schedule
  std::vector<Node> best_relinking_nodes;   //Best relinking active 
  bool one_improv = false;                  //Simple check to update just once the elite pool

  protected:

  /* assign_to_suboptimal:
  *   assign the given job to a suboptimal configuration, available in 
  *   an already opened node
  *
  *   Input:  const Job&                      job to be assigned
  *           const setup_time_t&             execution times of the given job
  *           Dstar&                          see dstar.hpp
  *           job_schedule_t&                 new proposed schedule
  *
  *   Output: bool                            true if job has been assigned
  */
  virtual bool assign_to_suboptimal (const Job&, const setup_time_t&,
                                     Dstar&, job_schedule_t&) override;

  /* perform_assignment 
  *   assign the selected job to the new proposed schedule
  *
  *   Input:  const Job&                      job to be assigned
  *           job_schedule_t&                 new proposed schedule
  */
  virtual void perform_assignment (const Job&, job_schedule_t&) override;

  /* perform_scheduling
  *   perform the scheduling process that assigns the best configuration to
  *   each submitted jobs (as long as resources are available)
  *
  *   Output: job_schedule_t        proposed schedule
  */
  virtual job_schedule_t perform_scheduling (void) override;

  /* scheduling_step
  *   sort the list of submitted jobs and perform scheduling of all
  *   submitted jobs
  *
  *   Input:    job_schedule_t&     empty schedule to be built
  */
  virtual void scheduling_step (job_schedule_t&) override;
  
  /* update_best_schedule by comparng the value of the proxy function
  *
  *   Input:    job_schedule_t&         new proposed schedule
  */
  void update_best_schedule (job_schedule_t&);

  /*  path relinking between two solutions:
  *    recursive function that explore the space
  *    until no moves between the current schedule
  *    and the guiding one or until the maximum number
  *    of moves is reached. At most one improvement is made
  *    in the elite pool (avoid homogeneization of the elite)
  *
  *   Input:    job_schedule_t&         source schedule
  *             const job_schedule_t&   guiding schedule
  *             double                  best improvement reached
  *             double                  current improvement             
  *             unsigned                depth of exploration
  * 
  *   Output:   Updated elite pool of solutions
  * 
  */
  void relinking_phase(job_schedule_t&, const job_schedule_t&, double, double, unsigned);

  /*  Explore the effect of applying a step to a schedule_t
  *    if the bool is set to true, the move is applied, else
  *    the schedule is reverted back to normal
  *    The functions also has a recursive part to search forbidden moves
  * 
  *   Input:    job_schedule_t&               Source schedule
  *             const job_schedule_t&         Target schedule       
  *             const pair<Job, int>&         Move to apply
  *             int                           Foresight (for recursion)
  *             bool                          Permament move (false by default)
  *
  *   Output:   double improvement of the modified schedule &
  *             modified schedule if the move is permanent
  * 
  */
  double explore_step(job_schedule_t &, const job_schedule_t&,
                      const std::pair<Job, int> &,
                      bool, bool permanent = false);

  /*  Search for possibile moves, given two schedules
  *    a move is valid if a node can accomodate
  *    the setup of the guiding solution, will
  *    generate a map containing all the possible moves
  *
  *   Input:  job_schedule_t& source
  *           job_schedule_t& guiding
  * 
  *   Output: unordered_map<Job, Int>
  * 
  */
  pool get_moves(job_schedule_t&, const job_schedule_t&);
  
  /*  check if two setups are identical, by GPU type 
  *   and number of GPU required
  *
  *   Input:  Setup& first setup
  *           Setup& second setup
  *     
  *   Output: bool
  * 
  */
  bool same_setup(const Schedule&, const Schedule&);

  /*  return a node index that is able to accomodate
  *    a proposed schedule. If no node is available
  *    return -1
  *   
  *   Input:  job_schedule_t& source schedule
  *           Job&            job to modify
  *           Schedule&       schedule to accomodate
  * 
  *   Output: int             node_index
  * 
  */
  int compatible(job_schedule_t &, const Job &, const Schedule &);
  
  /* This function just iterates trough the elite solutions
  *  and perform forward relinking between the best
  *  elite and another elite solution
  *  
  *  Output: modified best solution
  * 
  */
  void best_fight(void);

  /* Created to avoid computing each time the proxy function to test
  *   a move. It will apply the difference between the
  *   proxy before the move, and the proxy after the move.
  *   This results in computing the cost of the modified job,
  *   instead of the whole schedule. This function just a single
  *   iteration of greedy::obective_function()
  *
  * Input: const Job&       Job that will change setup
  *        const Schedule&  Outgoing schedule
  *        const Schedule&  Incoming schedule
  *        unsigned         GPUs of the node in the old Schedule      
  * 
  * Output: double          difference between the two proxy function
  * 
  */
  double update_best_cost(const Job&, const Schedule&, const Schedule&, unsigned);

  /* Print out some usefull infomation about 
  *  the path relinking method (in the output file)
  * 
  *  Input: double Greedy cost
  *         double RandomGreedy cost
  *         double PathRelinking cost 
  *   
  *  Output: current time, queue pressure,
  *          jobs in queue, nodes opened, 
  *          RG improvement, PR improvement
  *           
  */
  void print_info(double&, double&, double&);

  /*  print the moves in the output file (Job, VM and number of GPUs)
  *   
  *   Input : source and target schedule 
  *           pool of moves
  * 
  *   Output : list of moves in the output file
  */
  void print_moves(job_schedule_t &, const job_schedule_t &, pool &);

public:

  /*  constructor
  *
  *   Input:  const std::string&    ARGS = "nInitialJ-nN-nJ-lambdaa-mu-myseed"
  *           const std::string&    "-delta"
  *           const std::string&    name of file with the list of jobs
  *           const std::string&    name of file with execution times of jobs
  *           const std::string&    name of file with the list of nodes
  *
  */
  Path_relinking (const std::string&, const std::string&, 
                 const std::string&, const std::string&);
  virtual ~Path_relinking (void) = default;

};

