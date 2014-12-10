#include "caterpillar.hpp"

namespace status{

  unordered_map<cat_dynprog_input_pair, cat_dynprog_configs, input_pair_hasher> cat_DP_table;
  bool cat_DP_table_initialized(false);
  size_t largest_set_list(0);

  inline bool operator==(const cat_dynprog_input& X, const cat_dynprog_input& Y){
    return (X.status == Y.status) && (X.subtree == Y.subtree) && (X.influx == Y.influx);
  }
  inline bool operator==(const cat_dynprog_input_pair& X, const cat_dynprog_input_pair& Y){
    return (X.first == Y.first) && (X.second == Y.second);
  }

  inline vertex* add_new_leaf_with_status(tree* t, vertex* const parent, list<uint>& stati){
    vertex* result = t->add_vertex(parent);
    stati.pop_front();
    return result;
  }



  inline void mark_invalid(cat_dynprog_input& input){
    input = INVALID_INPUT;
  }
  inline bool is_invalid(const cat_dynprog_input& input){
    return input == INVALID_INPUT;
  }
  // lower bound on the influx based on number of vertices in the other side
  inline uint lower_influx_bound(const uint other_subtree){
    return other_subtree;
  }
  // upper bound on the influx based on number of vertices in the other side - when it forms a path - (thx to Euler ^^)
  inline uint upper_influx_bound(const uint other_subtree){
    return (other_subtree * (other_subtree + 1)) >> 1;
  }

  bool is_sane(const cat_dynprog_input& input, const sequence_t& seq, const uint num_vertices){
    // an input is insane if...
    DEBUG1(cout << "checking "<<input<<" for sanity..."<<endl);
    // ... not NO_INPUT
    if(input == NO_INPUT) return true;
    // ... its status is not in the sequence
    if(seq.find(input.status) == seq.end()) return false;
    DEBUG1(cout << input.status<<" passed containment in sequence"<<endl);
    // ... its subtree is larger than the graph
    if(input.subtree > num_vertices) return false;
    if(input.subtree == 0) return false;
    DEBUG1(cout << input.subtree<<" passed subtree bounds"<<endl);
    // ... its influx is less/more than the lower/upper bound based on vertices on the other side
    const uint other_side(num_vertices - input.subtree);
    if(input.influx < lower_influx_bound(other_side)) return false;
    if(input.influx > upper_influx_bound(other_side)) return false;
    DEBUG1(cout << input.influx<<" passed influx bounds..."<<endl);
    // ... its influx and subtree must generate a smaller/larger status
    if(input.influx + lower_influx_bound(input.subtree - 1) > input.status) return false;
    if(input.influx + upper_influx_bound(input.subtree - 1) < input.status) return false;
    DEBUG1(cout << input.status<< " passed status bounds..."<<endl);
    // if we're well within all limits, the input is sane
    return true;
  }

  // get next larger adjacent status using s' - s = |V| - 2n'
  inline uint get_adjacent_status(const uint my_status, const uint my_subtree_size, const uint num_vertices){
    return my_status + 2 * my_subtree_size - num_vertices;
  }

  // compute the next input on the path and add it to consumed_stati
  inline cat_dynprog_input compute_next_input(const cat_dynprog_input& old_input, const uint new_status, const uint leaves){
    cat_dynprog_input result;
    result.status = new_status;
    // compute the subtree using n = n' + l + 1
    result.subtree = old_input.subtree + leaves + 1;
    // compute the influx using g' + f' = s', g = g' + n' + l, and g + f = s together yielding f = f' + s - s' - n' - l
    result.influx = old_input.influx + result.status - old_input.status - old_input.subtree - leaves;
    return result;
  }

  // compute the infos of the new question given its leaves and the old inputs
  // if there are no leaves, then we assume we're looking for the center!
  cat_dynprog_input compute_inner_input(const cat_dynprog_input& old_input,
                                        const uint leaves,
                                        const uint sought_leaf_status,
                                        const uint num_vertices,
                                        list<uint>& consumed_stati){
    cat_dynprog_input result;
    // we expect to see the corresponding backbone status
    const uint expected_status( (leaves > 0) ? get_corresponding_backbone_status(sought_leaf_status, num_vertices) : sought_leaf_status);
    DEBUG2(cout << "walking inwards from "<<old_input<<", expecting to see "<<expected_status<<endl);
    // if there is NO_INPUT, then create a new input
    if(old_input == NO_INPUT){
      // create a new ending of the caterpillar
      return (cat_dynprog_input){
        expected_status,          // status
        leaves + 1,               // subtree
        expected_status - leaves  // influx
      };
    } else {
      // if we're looking for the backbone status that's outwards, just return the current input
      if((leaves == 0) && (expected_status >= old_input.status)) return old_input;

      // consume the vertex we're currently at
      consumed_stati.push_back(old_input.status);

      // compute the next status
      result.status = get_adjacent_status(old_input.status, old_input.subtree, num_vertices);
      DEBUG2(cout << "next status: " << result.status << ", expected "<<expected_status<<endl);

      // if the status didn't decrease, we went over the center
      if(result.status >= old_input.status) return INVALID_INPUT;
      // if the status decreased below what we're expecing, then the expected status cannot be reached
      if(result.status < expected_status) return INVALID_INPUT;

      // if the status matches the expected status, then we hit the jackpot :)
      if(result.status == expected_status) return compute_next_input(old_input, result.status, leaves); else {
        // if the status doesn't match what we're expecting, then we have to walk further
        result = compute_next_input(old_input, result.status, 0);
        return compute_inner_input(result, leaves, sought_leaf_status, num_vertices, consumed_stati);
      }
    }
  }

  
  
  // the inverse of the above
  cat_dynprog_input compute_outer_input(const cat_dynprog_input& old_input, const uint leaves, const uint num_vertices){
    cat_dynprog_input result;
    // compute the subtree using n = n' + l + 1
    result.subtree = old_input.subtree - leaves - 1;
    // compute the new status using s' - s = |V| - 2n'
    result.status = old_input.status - 2 * old_input.subtree + num_vertices;
    // compute the influx using g' + f' = s', g = g' + n' + l, and g + f = s together yielding f = f' + s - s' - n' - l
    result.influx = old_input.influx - result.status + old_input.status + old_input.subtree + leaves;
    return result;
  }


  // TODO: use (leaf stati <= exists vertex with n-2 status difference) to save guesses
  // initialize the dynamic programming table
  void caterpillar_dynprog_initialize(const sequence_t& stati_seq,
                                      const list<uint>& stati_list,
                                      const uint num_vertices)
  {
    DEBUG2(cout << "initializing dynprog table"<<endl);
    // initialize the set list counter as well
    largest_set_list = 0;
    // get the center status
    const uint center_status(stati_list.front());
    // get occurances of the center
    const uint center_occurances(stati_seq.at(center_status));
    // get leaf status corresponding to the center
    //const uint center_leaf_status(get_corresponding_leaf_status(center_status, num_vertices));
    // create container to hold the table entries
    cat_dynprog_config config;
    // leave stati_used blank

    // slightly hacky, please forgive me:
    // we're counting 0 each time we use the center, so we put the right number here
//    config.stati_used[center_status] = center_occurances;


    // the input stati to be computed
    cat_dynprog_input_pair inputs;
    if(center_occurances == 1){
      // the center is unique


      // for each partition of the vertices to the left and right side (by symmetry of the caterpillar, we assume left < right)
      for(uint subtree_left = 1; subtree_left <= (num_vertices + 1)/2; ++subtree_left){ 
        // if the left subtree is 0, then the second input is fully determined, go create the DP-table entry
        inputs.first.status = inputs.second.status = center_status;

        const uint subtree_right(num_vertices - subtree_left + 1);
        //const uint center_leaves((stati_seq.find(center_leaf_status) != stati_seq.end()) ? stati_seq.at(center_leaf_status) : 0);

        inputs.first.subtree = subtree_left;
        inputs.second.subtree = subtree_right;
        // influx to the right must be at least subtree_left-1
        DEBUG2(cout << "range for influx_right: "<<subtree_left-1<<"-"<<center_status-(subtree_right -1)<<endl);
        // for each partition of influxes and subgraph sizes, compute respective inputs (the root is counted twice)
        for(uint influx_right = subtree_left - 1; influx_right + (subtree_right - 1)<= center_status; ++influx_right){
          const uint influx_left(center_status - influx_right);
          inputs.second.influx = influx_right;
          inputs.first.influx = influx_left;

          // TODO: some of the inputs are still not sane (according to is_sane) - generate less!
          config.inner_cat = new tree();
//          vertex* root = config.inner_cat->add_vertex(NULL);
//          config.docks = make_pair(root, root);



          // finally, create a cat_DP_table entry with inputs and all possible status combinations on each side
          cat_dynprog_configs& confs(cat_DP_table[inputs]);
          confs.insert(config);
          DEBUG2(cout << "guess: subtree_left="<<subtree_left<<" of "<<num_vertices<<" influx_left="<<influx_left<<"\t\tT"<<inputs<<"="<<confs<<" now"<<endl);
        }
      }
    } else {
      // the center is not unique
      if(center_occurances > 2) FAIL("detected "<<center_occurances<<">2 center vertices");
      if(num_vertices & 1) FAIL("2 centers and odd number of vertices, this is not right...");

      // both sides get the same status and number of vertices by n - 2n_1 = s - s = n - 2n_2
      inputs.first.status = inputs.second.status = center_status;
      inputs.first.subtree = inputs.second.subtree = num_vertices>>1;

      // for each reasonable partition of influxes and leaves, compute respective inputs
      for(uint influx_left = 2*inputs.second.subtree - 1; influx_left + inputs.first.subtree - 1 <= center_status; ++influx_left){
        inputs.first.influx = influx_left;
        // use formulars f_2 + g_2 + l_2 = s  and  f_1 = g_2 +  l_2 + n_2
        inputs.second.influx = center_status + inputs.second.subtree - influx_left;

        config.inner_cat = new tree();
//        vertex* root1 = config.inner_cat->add_vertex(NULL);
//        vertex* root2 = config.inner_cat->add_vertex(root1);
//        config.docks = make_pair(root1, root2);

        // finally, create a cat_DP_table entry with inputs and all possible status combinations on each side
        cat_dynprog_configs& confs(cat_DP_table[inputs]);
        confs.insert(config);
        DEBUG2(cout << "guess: influx_left="<<influx_left<<"\t\tT"<<inputs<<"="<<confs<<" now"<<endl);
      }
    }
    // finally mark table as initialized
    cat_DP_table_initialized = true;
  }

  // update one half of a config given by a recursive call of the dynamic programming
  // return whether the forced stati are supported by the sequence
  void update_config_tree(tree*& t,
                          vertex*& dock,
                          const uint leaves)
  {
    // update the internal caterpillar t & the docking vertices if there are leaves
    DEBUG1(cout << "adding "<<leaves<<" leaves and a backbone vertex to tree @"<<t<<":"<<endl<< *t<<endl);
    // add the backbone vertex
    dock = t->add_vertex(dock);
    // add leaves to the last dock
    for(uint i = 0; i < leaves; ++i) t->add_vertex(dock);
    DEBUG1(cout << "yielded "<<endl<< *t<<endl);
  }

  // update the attachment with a given list of consumed stati
  bool update_attachments(sequence_t& stati_used,
                          const sequence_t& stati_seq,
                          const uint bb_status,
                          const uint leaf_status,
                          const uint center_status,
                          const uint leaves)
  {
    // add the backbone status to the attachment if we guessed so (otherwise it will be 0)
    // if we used more stati then we have, return failure
    if(leaves) {
      stati_used[leaf_status] += leaves;
      if(stati_seq.at(leaf_status) < stati_used.at(leaf_status)) return false;
    }

    stati_used[bb_status]++;
    if(stati_seq.at(bb_status) < stati_used.at(bb_status)) return false;

    return true;
  }

  // update a config given by a recursive call of the dynamic programming
  pair<bool, cat_dynprog_config> update_config(const cat_dynprog_config& c,
                     const sequence_t& stati_seq,
                     const pair<bool, bool>& update_who,
                     const pair<leaf_guess_t, leaf_guess_t>& guess,
                     const uint next_status,
                     const uint leaf_status,
                     const uint center_status)
  {
    // prepare container to hold the result
    pair<bool, cat_dynprog_config> result(true, c);
    cat_dynprog_config& rc(result.second);

    if(update_who.first){
      if((next_status == center_status) && (stati_seq.at(center_status) == 1) && (guess.first.leaves == 0)){
        DEBUG2(cout << "skipping update for virtual left center"<<endl);
      } else {
        if(!update_attachments(rc.stati_used, stati_seq, next_status, leaf_status, center_status, guess.first.leaves)){
          result.first = false;
          return result;
        }
        update_config_tree(rc.inner_cat, rc.docks.first, guess.first.leaves);
      }
    }
    if(update_who.second){
      if(!update_attachments(rc.stati_used, stati_seq, next_status, leaf_status, center_status, guess.second.leaves)){
        result.first = false;
        return result;
      }
      update_config_tree(rc.inner_cat, rc.docks.second, guess.second.leaves);
    }
    return result;
  }


  // forward-declare dynamic programming
  cat_dynprog_configs caterpillar_dynprog(const sequence_t&,
                                          const list<uint>&,
                                          const cat_dynprog_input_pair&,
                                          const uint);

  cat_dynprog_configs dynprog_recurse_for_guess(const sequence_t& stati_seq,
                                                const list<uint>& stati_list,
                                                const uint next_status,
                                                const cat_dynprog_input_pair& inputs,
                                                pair<bool, bool> update_who,
                                                const pair<leaf_guess_t, leaf_guess_t>& guess,
                                                const uint num_vertices)
  {
    cat_dynprog_input_pair new_inputs(inputs);
    
    // 3. take 1 step inwards to get the new input values
    const uint center_status(stati_list.front());
    const uint leaf_status(get_corresponding_leaf_status(next_status, num_vertices));

    // don't update the first input if it is NO_INPUT and the guess says that it's not getting any leaves
    if((inputs.first == NO_INPUT) && (guess.first.leaves == 0)) update_who.first = false;

    if(update_who.first){
      if(inputs.first == NO_INPUT)
        new_inputs.first = (cat_dynprog_input){
                      next_status,                     // status
                      guess.first.leaves + 1,          // subtree
                      next_status - guess.first.leaves // influx
                    };
      else {
        new_inputs.first  = compute_next_input(inputs.first,  next_status, guess.first.leaves);
      }
    }
    if(update_who.second){
      new_inputs.second = compute_next_input(inputs.second, next_status, guess.second.leaves);
    }

    DEBUG2(cout<<"updating "<<update_who<<" yielded new inputs="<<new_inputs<<" from "<<inputs<<endl);

    // prepare container to hold the result
    cat_dynprog_configs result;
    // make sure the inputs are sane
    if(is_sane(new_inputs.first, stati_seq, num_vertices) && is_sane(new_inputs.second, stati_seq, num_vertices)){
      // recursive call to the previous layer of the dynamic programming table
      cat_dynprog_configs tmp(caterpillar_dynprog(stati_seq, stati_list, new_inputs, num_vertices));

      // remove all configurations that do not support our old inputs
      for(cat_dynprog_config c : tmp){
        DEBUG2(cout << "used so far: "<<c.stati_used<<", want to add: "<<next_status<<" & "<<guess.first.leaves+guess.second.leaves<<'x'<<leaf_status<<endl);
        pair<bool, cat_dynprog_config> updated(update_config(c, stati_seq, update_who, guess, next_status, leaf_status, center_status));
        if(updated.first){
          DEBUG2(cout << "supported, now used "<< updated.second.stati_used<<endl);
          // if we're at the top level (indicated by input NO_INPUT,NO_INPUT), no more stati should be attached
          if((inputs.first == NO_INPUT) && (inputs.second == NO_INPUT)){
            if(equal(stati_seq, updated.second.stati_used)) result.insert(updated.second); else
              DEBUG2(cout<< updated.second.stati_used << " does not exactly use our stati, discarding"<<endl);
          } else result.insert(updated.second);
        } else DEBUG2(cout << "not supported!"<<endl);
          // if the forced stati don't match the lists supported by the dynamic table entries, then erase this set-list
      }
      return result;
    } else return cat_dynprog_configs();
  }

  // do the dynamic programming
  cat_dynprog_configs caterpillar_dynprog(const sequence_t& stati_seq,
                                          const list<uint>& stati_list,
                                          const cat_dynprog_input_pair& inputs,
                                          const uint num_vertices)
  {
    // if cat_DP_table has not yet been initialized, do so first
    if(!cat_DP_table_initialized) caterpillar_dynprog_initialize(stati_seq, stati_list, num_vertices);
    // if it's in the DP table, use it
    { const auto lookup(cat_DP_table.find(inputs));
      if(lookup != cat_DP_table.end()){
        DEBUG2(cout << "=== found "<<inputs<<" with "<<lookup->second.size()<<" entries in the table:"<<endl);
        return lookup->second;
    }}
    DEBUG2(cout << inputs << " not found in the table, computing..."<<endl);
    // get the center status
    const uint center_status(stati_list.front());
    // prepare container to hold result
    cat_dynprog_configs result;
    // the next status
    uint next_status;
    // who is going to be included in the guesswork
    pair<bool, bool> advance_who(make_pair(false, true));

    // compute the next status and who is going to be advanced
    if(inputs.second == NO_INPUT){
      next_status = get_corresponding_backbone_status(stati_list.back(), num_vertices);
      advance_who.first = true;
    } else {
      next_status = get_adjacent_status(inputs.second.status, inputs.second.subtree, num_vertices);
      DEBUG2(cout<<" computed "<<next_status<<" from "<<inputs.second.status<<" and "<<inputs.second.subtree<<endl);
      if(inputs.first == NO_INPUT){
        if(inputs.second.status == center_status){
          // if we're at the center, but inputs.first == NO_INPUT, then just give it [x 1 x], where x is the canter status
          result = caterpillar_dynprog(stati_seq, stati_list,
              make_pair( (cat_dynprog_input){
                  center_status, // status
                  1,             // subtree
                  center_status  // influx
                }, inputs.second),
              num_vertices);
          // register the center as used in all configs
          for(auto c : result) c.stati_used[center_status]++;
          return result;
        } else {
          // if we found something invalid, return failure
          if((next_status > inputs.second.status) || (stati_seq.find(next_status) == stati_seq.end())) return cat_dynprog_configs();
          // if the first input is NO_INPUT but the second is not (and its not the center)
          // include the first input in the guesswork (if the sequence supports another backbone vertex of this status)
          if(stati_seq.at(next_status) > 1) advance_who.first = true;
          // recurse for every input status strictly between the first and its adjacent next one
          DEBUG2(cout << "checking if anyone between "<<next_status<<" and "<<inputs.second.status<<" in "<<stati_list<<" could be at the left end"<<endl);
          for(list<uint>::const_iterator i = stati_list.begin(); *i < inputs.second.status; ++i) if(*i > next_status) {
            DEBUG2(cout << *i << " is a candidate for the left end"<<endl);
            const uint leaf_status(get_corresponding_leaf_status(*i, num_vertices));
            const sequence_t::const_iterator leaf_occurances = stati_seq.find(leaf_status);
            // if there is no corresponding leaf-status, then this status is not what we're looking for (both ending should have leaves)
            if(leaf_occurances == stati_seq.end()) continue;

            DEBUG2(cout << "since we have only one input, I'll recurse for "<<*i<<" with leaf status "<<leaf_status<<" occuring "<<leaf_occurances->second<<"x"<<endl);
            if(leaf_occurances != stati_seq.end()){ // construct recursive calls for 0 or 1 occurance on the right backbone
              unordered_set_union(result, dynprog_recurse_for_guess(stati_seq, stati_list, *i, inputs, make_pair(true, false),
                    make_pair( (leaf_guess_t){ 0, leaf_occurances->second }, NEW_GUESS), num_vertices));
              // maybe leaf_status occurs also on the backbone on the right, so try using 1 leaf less
              if(leaf_occurances->second > 1)
                unordered_set_union(result, dynprog_recurse_for_guess(stati_seq, stati_list, *i, inputs, make_pair(true, false),
                      make_pair( (leaf_guess_t){ 0, leaf_occurances->second - 1 }, (leaf_guess_t){ 1, 0 }), num_vertices));
            }

          }
        }
      } else {
        // if none of the input stati is NO_INPUT, but we're over the center, then return failure
        if(next_status >= inputs.second.status) return cat_dynprog_configs();
        // also, if the status we computed doesn't exist, return failure
        if(stati_seq.find(next_status) == stati_seq.end()) return cat_dynprog_configs();


        // check if the adjacent status on the left is larger than on the right
        const uint tmp(get_adjacent_status(inputs.first.status, inputs.first.subtree, num_vertices));
        // and if so, update next_status
        if(tmp > next_status){
          next_status = tmp;
          advance_who = make_pair(true, false);
        } else advance_who.first  = (tmp == next_status);
      }
    }
    
    // get the next leaf status
    const uint next_leaf_status(get_corresponding_leaf_status(next_status, num_vertices));
    // get the occurances of the next leaf status in the sequence
    const uint next_occurances( (stati_seq.find(next_leaf_status) != stati_seq.end() ? stati_seq.at(next_leaf_status) : 0));

    // 2. guess distribution of occurances of next_leaf_status
    // branch into partitions of the vertices of next status:
    // backbone left? backbone right? #leaves left? #leaves right?
    pair<leaf_guess_t, leaf_guess_t> guess;

    DEBUG2(cout << next_status << " corresponds to "<<next_leaf_status<<" which occurs "<<next_occurances<<"x"<<endl);
    const uint max_bb_first(min(1U, next_occurances));
    DEBUG2(cout << "max guess for first backbone: "<<max_bb_first<<endl);
    for(guess.first.backbone = 0; guess.first.backbone <= max_bb_first; ++guess.first.backbone){
      const uint max_bb_second(min(1U, next_occurances - guess.first.backbone));
      // note: if next_leaf_status is larger than both input stati, then we only have to try putting it on the "outside backbone" once
      unsigned char min_bb_second = 0;
      if(next_leaf_status > max(inputs.first.status, inputs.second.status))
        if(guess.first.backbone == 1) min_bb_second = 1;
      DEBUG2(cout << "max guess for second backbone: "<<max_bb_second<<endl);
      for(guess.second.backbone = min_bb_second; guess.second.backbone <= max_bb_second; ++guess.second.backbone){
        const uint num_leaves(next_occurances - guess.first.backbone - guess.second.backbone);
        DEBUG2(cout << "max guess for first leaves: "<< (advance_who.first ? num_leaves : 0) <<endl);
        for(guess.first.leaves = (advance_who.second ? 0 : num_leaves); guess.first.leaves <= ( advance_who.first ? num_leaves : 0); ++guess.first.leaves){
          guess.second.leaves = num_leaves - guess.first.leaves;
          DEBUG2(cout << "guessed distribution "<<guess<<" of "<<next_leaf_status<<endl);

          unordered_set_union(result, dynprog_recurse_for_guess(stati_seq, stati_list, next_status, inputs, advance_who, guess, num_vertices));
          DEBUG2(cout << "results for inputs "<<inputs<<" augmented to "<<result<<endl);
        }
      }
    }

    // add the result to the dynamic programming table
    cat_DP_table[inputs] = result;
    // keep track of list sizes
    largest_set_list = max(largest_set_list, result.size());

    return result;
  }

  tree* stati_to_caterpillar(const sequence_t& s){
    // we rather work with a (sorted) list of stati
    list<uint> stati(get_occuring_stati(s));
    stati.sort();

    DEBUG2(cout << "stati: "<<stati<<endl);

    // let's say we have at least two stati
    assert(stati.size() > 1);

    // construct the inputs
    cat_dynprog_input_pair inputs;
    pair<uint, uint> leaves;

    const uint num_vertices(get_num_vertices(s));

    cat_dynprog_configs confs(caterpillar_dynprog(s, stati, cat_dynprog_input_pair(NO_INPUT, NO_INPUT), num_vertices));
    // and, if successfull, return the first possible tree
    if(!confs.empty()) return confs.begin()->inner_cat; else return NULL;
  }

  size_t get_set_list_max(){ return largest_set_list; }
}

std::ostream& operator<<(std::ostream& os, const status::cat_dynprog_input& i){
  return os << '['<<i.status<<' '<<i.subtree<<' '<<i.influx<<']';
}


ostream& operator<<(ostream& os, const status::cat_dynprog_config& c){
  os << c.stati_used;
  DEBUG1(os<< " tree @"<<c.inner_cat);
  return os;
}


ostream& operator<<(ostream& os, const status::leaf_guess_t& g){
  return os << '('<<(uint)g.backbone<<','<<g.leaves<<')';
}
