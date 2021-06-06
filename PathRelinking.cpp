#include "path_relinking.hpp"

Path_relinking::Path_relinking (const std::string& directory, 
                              const std::string& file_jobs, 
                              const std::string& file_times, 
                              const std::string& file_nodes):
  Random_greedy(directory, file_jobs, file_times, file_nodes)
{}

//The following functions are from the Random Greedy class, only
//minor (nontheless important) modifications have been made to them

job_schedule_t
Path_relinking::perform_scheduling(){
    
  //The pool of elite solution is empty at first
  K_best.clear();
  //We start by performing the random greedy procedure:
  preprocessing();
  // initialization of minimum total cost, best schedule and corresponding
  // index by a step of pure greedy
  job_schedule_t full_greedy_schedule;
  Greedy::scheduling_step(full_greedy_schedule);
  update_best_schedule(full_greedy_schedule);
  //print full greedy cost
  double cost_FG = K_best.rbegin()->first;
  std::cout << ", "<< cost_FG;
  
  //Random Greedy scheduling
  for (unsigned random_iter = 1; random_iter < (max_random_iter); 
       ++random_iter)
  {
    //generate the new schedule
    job_schedule_t new_schedule;
    scheduling_step(new_schedule);
    //Update the elite pool if possibile 
    update_best_schedule(new_schedule);
  }

  //print the cost after the random greedy
  double cost_RG = K_best.rbegin()->first;
  std::cout << ", " << cost_RG;

  //CORE OF THE PR: perform the PR between elite solutions
  //This function can be moved or put in a cycle to perform
  //different types of path relinking.
  best_fight();

  //print the cost after the path relinking
  double cost_PR = K_best.rbegin()->first;
  std::cout << ", " << cost_PR;

  //get best solution
  Best_solution& BS = K_best.rbegin()->second;
  job_schedule_t BS_sch = BS.get_schedule();
  std::swap(nodes,BS.get_opened_nodes());
  last_node_idx = BS.get_last_node_idx();

  //Print infos about the time step
  print_info(cost_FG, cost_RG, cost_PR);

  //Return the best solution
  return BS_sch;
};

void
Path_relinking::scheduling_step (job_schedule_t& new_schedule)
{
  if (r_swap)
  {
    sort_jobs_list();
    random_swap();
  }
  
  close_nodes();

  std::string queue = "";

  // perform assignment of all submitted jobs (STEP #1)
  // Note, I've added a check to see if all the resources are
  // taken or not (save a lot of time in exponential scenario)
  bool resources = true;
  for (const Job&j : submitted_jobs)
  {
    queue += (j.get_ID() + "; ");
    if(resources){
      perform_assignment(j, new_schedule);
      resources = check_resources();
    }
    else{
      new_schedule[j]=Schedule();
    }
  }
  // perform postprocessing (STEP #2)
  postprocessing(new_schedule);
}

void
Path_relinking::update_best_schedule (job_schedule_t& new_schedule)
{
  // find execution time of first ending job (historic reasons)
  double first_finish_time = find_first_finish_time(new_schedule);

  // compute cost of current schedule
  double current_cost = objective_function(new_schedule, first_finish_time,
                                           nodes);  

  // check if current solution is one of the best schedules
  K_best_t::iterator ub = K_best.upper_bound(current_cost);

  Best_solution BS(new_schedule, nodes, last_node_idx, first_finish_time);
  K_best.insert(ub, std::pair<double, Best_solution>(current_cost, BS));
  if (K_best.size() > K)
  {
    K_best_t::iterator it = K_best.begin();
    K_best.erase(it);
  }
}

bool
Path_relinking::assign_to_suboptimal (const Job& j, const setup_time_t& tjvg,
                                     Dstar& dstar, 
                                     job_schedule_t& new_schedule)
{
  bool assigned = false;
  while (!dstar.is_end() && !assigned)
  {
    dstar.set_generator(generator);
    setup_time_t::const_iterator best_stp = dstar.get_best_setup();
    generator = dstar.get_generator();
    assigned = assign_to_existing_node(j, best_stp, new_schedule);
  }
  return assigned;
}

void
Path_relinking::perform_assignment (const Job& j, job_schedule_t& new_schedule)
{
  // determine the setups s.t. deadline is matched
  const setup_time_t& tjvg = ttime.at(j.get_ID());
  Dstar dstar(j, tjvg, current_time);
  dstar.set_random_parameter(alpha);

  // determine the best setup:
  //   if it is possible to match deadline, the best setup is the 
  //   cheapest; otherwise, the best setup is the fastest
  setup_time_t::const_iterator best_stp;
  dstar.set_generator(generator);
  best_stp = dstar.get_best_setup();
  generator = dstar.get_generator();

  // check the already opened nodes...
  bool assigned = assign_to_existing_node(j, best_stp, new_schedule);

  // if the already opened nodes are not suitable...
  if (! assigned)
  {
    // if it is possible, open a new node with the optimal configuration
    // else, assign to an available suboptimal configuration
    if (last_node_idx < nodes.size())
    { 
      assigned = true;
      assign_to_new_node(j, best_stp, new_schedule);
    }
    else
      assigned = assign_to_suboptimal(j, tjvg, dstar, new_schedule);
  }

  // if job j cannot be assigned to any configuration, add an empty 
  // schedule
  if (!assigned)
    new_schedule[j] = Schedule(); 
}

//From here, there are all functions are created for the project.

void
Path_relinking::best_fight(){

  K_best_t Kb = K_best; //We copy the starting elite pool

  for (auto it = Kb.begin(); it != Kb.end(); it++)
  {
    Best_solution &FB = K_best.rbegin()->second;            //Best solution
    nodes = FB.get_opened_nodes();                          //Best schedule nodes
    job_schedule_t FB_sch = FB.get_schedule();              //Best schedule
    job_schedule_t &next_elite = it->second.get_schedule(); //Next elite

    //Path relinking
    explored_moves.clear();
    relinking_phase(FB_sch, next_elite, 0., 0., MAX_DEPTH); //Perform the relinking
  }
}

void
Path_relinking::relinking_phase(job_schedule_t &source, 
                                const job_schedule_t &target, 
                                double relinking_best, double improvement,
                                unsigned DEPTH)
{
  std::pair<double, std::pair<Job, int>> step_best = {}; //Best improv & move
  pool moves = get_moves(source, target);                //Get legal moves     
  double fitness=0;                                      //Initial fitness value       

  //If a move is possible, perform the linking
  if((!moves.empty()) && DEPTH>0){
    //Consume a relinking move
    DEPTH--; 
    //Explore the moves (not permanently)
    for (auto move : moves)
    {          
      //Check if the move has already been explored
      auto fit = explored_moves.find(move);
      if(fit==explored_moves.end()){
        fitness = explore_step(source, target, move, FUTURE_SIGHT);
        explored_moves[move] = fitness;
      }else{
        fitness = fit->second;
      }
      //Find the best step fitness & the associated move!
      if (fitness >= step_best.first || step_best.first==0){
        step_best = {fitness, move};
      }
    }
    //Apply the best move (Permanently)
    improvement += explore_step(source, target, step_best.second, false, true);  

    //Update the best solution found in the linking if possible
    if (improvement >= relinking_best)
    {
      best_relinking_sch = source;       //Schedule of the best solution
      best_relinking_nodes = nodes;      //Nodes of the best solution
      relinking_best = improvement;      //Update the best relinking found
      one_improv = true;                 //One improvement is found!
    }
    //Move to the next step! Recursion!
    relinking_phase(source, target, relinking_best, improvement, DEPTH); 
  }
  //If an improvement has been found during the relinking, update 
  if(one_improv)
  {
    one_improv=false;                         //Improve only once
    std::swap(nodes,best_relinking_nodes);    //Select the best nodes
    postprocessing(best_relinking_sch);       //Apply post-processing
    update_best_schedule(best_relinking_sch); //Update the elite pool
  }
}

double
Path_relinking::explore_step(job_schedule_t &source, const job_schedule_t &target, 
                              const std::pair<Job, int> &move, 
                              bool FORESIGHT, bool permanent){

  Job j = move.first;                //Job to be modified
  Schedule ins = target.at(j);       //Setup to insert
  Schedule rem = source.at(j);       //Setup to remove

  unsigned node_idx = 0;             //Initialization of node index  
  unsigned GPU_source = 0;           //Initialization of GPU required

  unsigned new_node_idx = move.second;   //Node of the new setup
  int GPU_target = 0;                    //GPUs of the new setup
  unsigned GPU_rem = 0;                  //Total GPUs in the node

  //Remove the old setup
  if (!rem.isEmpty())
  {
    GPU_source = source[j].get_setup().get_nGPUs();     //GPUs to free
    node_idx = source[j].get_node_idx();                //Node containing the job
    GPU_rem = nodes[node_idx].get_usedGPUs();           //GPUs used before removal
    nodes[node_idx].set_remainingGPUs(-GPU_source);     //GPUs are freed from the node

    //If there are no more jobs in the node, we close it
    if (nodes[node_idx].get_usedGPUs() == 0){
      nodes[node_idx].close_node();
    }
  }

  //Insert the new setup
  if (!ins.isEmpty())
  {
    GPU_target = ins.get_setup().get_nGPUs();          //GPUs of the new setup
    if (!nodes[new_node_idx].open()){                  //Open the node if necessary
      nodes[new_node_idx].open_node(ins.get_setup());
    }
    nodes[new_node_idx].set_remainingGPUs(GPU_target); //GPUs are activated
    ins.set_node(new_node_idx);                        //Setting the node index
  }

  //Modify the Schedule
  source[j] = ins; 

  //Check what is the difference in old and new proxy function 
  double cost = update_best_cost(j,ins,rem, GPU_rem); 

  if(FORESIGHT){
    //Calculate the new moves
    pool moves_FS = get_moves(source, target); 
    //Remove moves that were already available before
    for (auto it = moves.begin(); it != moves.end(); it++){
      moves_FS.remove(*it);
    }
    //The base has to stay the same each iteration
    double baseline = cost; 
    for (auto move : moves_FS){
      //Explore one step ahead
      double price = baseline + explore_step(source, target, move, false);
      if (price > cost){
        cost = price; //Update the cost of the move with the best one
      }
    }
  }

  //Revert the situation to the previous one if not permanent
  if (!permanent){
    //Guiding setup is removes
    if(!ins.isEmpty()){
      nodes[new_node_idx].set_remainingGPUs(-GPU_target); //GPUs are freed from the node
      if (nodes[new_node_idx].get_usedGPUs() == 0){       //Check the only job condition
        nodes[new_node_idx].close_node();                 //Node is resetted
      }
    }
    
    //Starting setup is inserted
    if (!rem.isEmpty()){
      if (!nodes[node_idx].open()){                      //Open the node if necessary
        nodes[node_idx].open_node(rem.get_setup());      //Configure the node
      }
      nodes[node_idx].set_remainingGPUs(GPU_source);     //GPUs are activated
    }
    //Schedule is reverted back to normal
    source[j] = rem; 
  }

  //Cost is returned
  return cost; 
}

pool 
Path_relinking::get_moves(job_schedule_t &source, const job_schedule_t &target)
{
  pool c_moves = {};  //Pool of moves
  int node_idx = {};  //Node idx

  for(Job& j: submitted_jobs){
    if (!target.at(j).isEmpty())
    {
      //Setups must be different
      if (!same_setup(source.at(j), target.at(j))){
        //Search for a node that can accomodate the setup
        node_idx = compatible(source, j, target.at(j)); 
        if (node_idx!=-1){
          //Add move to the pool
          c_moves.push_back({j, node_idx});
        };
      };
    }
    //There is always space to remove a Setup
    else if(!source.at(j).isEmpty() && target.at(j).isEmpty()){
      c_moves.push_back({j, -1}); //-1 is just a placeholder
    }
  };
  //Print the moves in the output file (for deubbing)
  //print_moves(source, target, moves);
  return c_moves;
};

int
Path_relinking::compatible(job_schedule_t &source, const Job &j, const Schedule &candidate_schedule){

  unsigned Required_GPU = candidate_schedule.get_setup().get_nGPUs(); //GPUs needed by the new setup
  unsigned node_idx = 0;                                              //Index of the new setup

  //First check the node the old setup was at:
  if (!source[j].isEmpty()){  
    node_idx = source[j].get_node_idx();
    //If the job was the only one in the node, we can accomodate the new setup for sure
    if (nodes[node_idx].get_usedGPUs() == source[j].get_setup().get_nGPUs()){
      return node_idx;
    }
    //If the job is in a shared node, we search if there is space for the new setup
    //Note: now the VM type is important here as it cannot be changed
    else{
      //Check VM and GPU type
      if (source[j].get_setup() == candidate_schedule.get_setup()){
        //Calculate if enough GPUs are available for the new setup
        int GPU_avail = nodes[node_idx].get_remainingGPUs();                           
        int GPU_diff  = source[j].get_setup().get_nGPUs() + GPU_avail - Required_GPU;  
        if (GPU_diff >= 0){
          return node_idx;
        }
      }
    }
  }

  //If the old node is not available, search in the remaining ones
  node_idx = 0;
  for (Node& node : nodes){
    if (node.open()){
      //3 checks: Same VM, Same GPU type and enough GPU are free to accomodate the new setup
      bool same_VM    = node.get_VMtype() == candidate_schedule.get_setup().get_VMtype();   
      bool same_GPU   = node.get_GPUtype() == candidate_schedule.get_setup().get_GPUtype(); 
      bool enough_GPU = node.get_remainingGPUs() >= Required_GPU;           

      if (same_VM && same_GPU && enough_GPU){
        //Return the index of the node that can host the setup
        return node_idx;
      }
    }
    else{
      //Return the idx of the closed node
      return node_idx;
    }
    //Check next node
    node_idx++; 
  }
      
  return -1; //Special value if no nodes have been found
};


bool
Path_relinking::same_setup(const Schedule &source_sch, const Schedule &target_sch){

  //There was already an equality operator for setups
  //that only compared the VM and GPU types
  if(!source_sch.isEmpty() && !target_sch.isEmpty()){
    if(source_sch.get_setup().get_nGPUs()==target_sch.get_setup().get_nGPUs()){
      return(source_sch.get_setup()==target_sch.get_setup());
    }
  }
  //If both schedules are empty, the "setup" is still considered the same
  else if(source_sch.isEmpty() && target_sch.isEmpty()){
    return true;
  }
  return false;
};

double
Path_relinking::update_best_cost(const Job& j, const Schedule& ins, 
                                 const Schedule& rem, unsigned used_GPU_rem){

  //This function is extrimely ugly, it just perform the difference
  //between two proxy functions, these are just operations, nothing
  //much interesting happens here

  double price_gain =0;

  if(!rem.isEmpty()){
    //Calculate the old setup price
    double finish_time = rem.get_selectedTime();
    unsigned GPUs_rem = rem.get_setup().get_nGPUs();
    double vm_cost = finish_time * rem.get_setup().get_cost() / 3600 * GPUs_rem/used_GPU_rem;
    double tardi_cost = j.get_tardinessWeight() * std::max(current_time + finish_time - j.get_deadline(), 0.0);

    price_gain = -j.get_maxExecTime() / (vm_cost+tardi_cost);
  }
  if (!ins.isEmpty()){
    //Calculate the new setup price
    unsigned node_idx = ins.get_node_idx();
    unsigned g = nodes[node_idx].get_usedGPUs();
    double finish_time = ins.get_selectedTime();
    unsigned GPUs_ins = ins.get_setup().get_nGPUs();
    double vm_cost = finish_time * ins.get_setup().get_cost() / 3600 * GPUs_ins / g;
    double tardi_cost = j.get_tardinessWeight() * std::max(current_time + finish_time - j.get_deadline(), 0.0);

    price_gain += j.get_maxExecTime() / (vm_cost + tardi_cost);
  }

  return price_gain;
}

//This functions are dedicated to print information in the output file ("build/output/")

void
Path_relinking::print_info(double& cost_FG, double& cost_RG, double& cost_PR){ 

  //Calculate the number of opened nodes
  unsigned opened_nodes=0;
  for(auto &node:nodes){
    //I just wanted to use a regex
    (node.open()) ? opened_nodes++ : 0;
  }

  //Calculate the job pressure, normalized by the
  //system size
  double job_p=0;
  for(auto &j: submitted_jobs){
    job_p += j.get_pressure()/nodes.size();
  }

  //Get the improvement from using PR:
  double improv = ((cost_RG+1)<cost_PR) ? cost_PR-cost_RG : 0;
  double improv_RG = ((cost_FG+1)<cost_RG) ? cost_RG-cost_FG : 0;

  //Print in the output file some usefull informations
  std::cout 
  << ", " << current_time           //Sim time
  << ", " << job_p                  //Job pressure
  << ", " << submitted_jobs.size()  //Queue size
  << ", " << opened_nodes           //Opened nodes
  << ", " << improv_RG              //RG improvement
  << ", " << improv                 //PR improvement
  << std::endl;
}

void 
Path_relinking::print_moves(job_schedule_t &source, const job_schedule_t &target, pool &moves)
{
  //This function just print the moves, usefull to check runtime errors during testing
  if (!moves.empty())
  {
    std::cout << std::endl;
    for (auto move : moves)
    {
      if (!source[move.first].isEmpty() && !target.at(move.first).isEmpty())
      {
        std::cout << move.first.get_ID() << " from: " << source[move.first].get_setup().get_VMtype()
                  << ", " << source[move.first].get_setup().get_nGPUs() << " to: " << target.at(move.first).get_setup().get_VMtype()
                  << ", " << target.at(move.first).get_setup().get_nGPUs() << std::endl;
      }
      else if(!target.at(move.first).isEmpty()){
        std::cout << move.first.get_ID() << " from: EMPTY SCHEDULE" << " to: " << target.at(move.first).get_setup().get_VMtype()
                  << ", " << target.at(move.first).get_setup().get_nGPUs() << std::endl;
      }
      else{
        std::cout << move.first.get_ID() << " from: " << source[move.first].get_setup().get_VMtype()
                  << ", " << source[move.first].get_setup().get_nGPUs() << " to: " << "EMPTY SCHEDULE" << std::endl;
      }
    }
    std::cout << std::endl;
  }
};
