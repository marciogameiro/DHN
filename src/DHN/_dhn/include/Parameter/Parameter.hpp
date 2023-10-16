/// Parameter.hpp
/// Shaun Harker
/// 2015-05-24
///
/// Marcio Gameiro
/// 2021-07-13

#pragma once

#ifndef INLINE_IF_HEADER_ONLY
#define INLINE_IF_HEADER_ONLY
#endif

#include "Parameter.h"

INLINE_IF_HEADER_ONLY Parameter::
Parameter ( void ) {
  data_ . reset ( new Parameter_ );
}

INLINE_IF_HEADER_ONLY Parameter::
Parameter ( std::vector<LogicParameter> const& logic,
            std::vector<OrderParameter> const& order,
            Network const& network ) {
  assign ( logic, order, network );
}

INLINE_IF_HEADER_ONLY Parameter::
Parameter ( Network const& network ) {
  assign ( network );
}

INLINE_IF_HEADER_ONLY void Parameter::
assign ( std::vector<LogicParameter> const& logic,
         std::vector<OrderParameter> const& order,
         Network const& network ) {
  data_ . reset ( new Parameter_ );
  data_ -> logic_ = logic;
  data_ -> order_ = order;
  data_ -> network_ = network;
}

INLINE_IF_HEADER_ONLY void Parameter::
assign ( Network const& network ) {
  data_ . reset ( new Parameter_ );
  data_ -> network_ = network;
}

INLINE_IF_HEADER_ONLY bool Parameter::
attracting ( Domain const& dom ) const {
  int D = data_ -> network_ . size ();
  for ( int d = 0; d < D; ++ d ) {
    if ( not dom . isMin(d) && absorbing ( dom, d, -1 ) ) return false;
    if ( not dom . isMax(d) && absorbing ( dom, d, 1 ) ) return false;
  }
  return true;
}

INLINE_IF_HEADER_ONLY std::vector<bool> Parameter::
combination ( Domain const& dom, int variable ) const {
  std::vector<bool> input_combination;
  //std::cout << "  Forming input combination by analyzing inputs of node " << variable << ".\n";
  for ( int source : data_ -> network_ . inputs ( variable ) ) {
    //std::cout << "    Analyze source edge " << source << "\n";
    bool activating = data_ -> network_ . interaction ( source, variable );
    //std::cout << "      This edge is " << (activating ? "activating" : "repressing" ) << ".\n";
    int inedge = data_ -> network_ . order ( source, variable );
    //std::cout << "      This edge is the " << inedge << "th ordered outedge of " << source << ".\n";
    int thres = data_ -> order_ [ source ] . inverse ( inedge );
    //std::cout << "      The input combination digit depends on which side of threshold " << thres << " on dimension " << source << " we are at.\n";
    bool result = not ( dom [ source ] > thres ) ^ activating;
    //std::cout << "      The domain is on the " << ( ( dom [ source ] > thres ) ? "right" : "left" ) << " side of this threshold.\n";
    input_combination . push_back ( result );
    //std::cout << "      Hence, the input combination digit is " << (result ? "1" : "0") << "\n";
  } 
  //std::cout << "  Input combination formed. Big-endian representation = ";
  //for ( int i = input_combination . size () - 1; i >= 0; -- i ){
  // std::cout << (input_combination[i] ? "1" : "0");
  //}
  //std::cout << "\n";
  return input_combination;
}

INLINE_IF_HEADER_ONLY bool Parameter::
absorbing ( Domain const& dom, int collapse_dim, int direction ) const {
  //std::cout << "Absorbing (" << dom . index () << ", " << collapse_dim << ", " << direction << ")\n";
  int thres = dom [ collapse_dim ];
  if ( direction == -1 ) thres -= 1;
  if ( thres < 0 ) return false;
  if ( thres == data_ -> network_ . num_thresholds (collapse_dim) ) return false;
  // if ( thres == data_ -> network_ . outputs( collapse_dim).size() ) return false;
  //std::cout << "  Threshold # = " << thres << "\n";
  std::vector<bool> input_combination = combination(dom, collapse_dim);
  //std::cout << "  Consulting parameter " <<  data_ -> logic_ [ collapse_dim ] . stringify () << ".\n";
  bool flow_direction = data_ -> logic_ [ collapse_dim ] ( input_combination, thres );
  //std::cout << "  Flow direction is to the " << (flow_direction ? "right" : "left") << "\n";
  if ( direction == -1 ) {
    //std::cout << "  Hence the left wall is " << ((flow_direction)?"not ":"") << "absorbing.\n";
    return not flow_direction;
  } else {
    //std::cout << "  Hence the right wall is " << ((not flow_direction)?"not ":"") << "absorbing.\n";
    return flow_direction;
  }
}

INLINE_IF_HEADER_ONLY uint64_t Parameter::
regulator ( uint64_t variable, uint64_t threshold ) const {
  uint64_t num_outedges = network () . outputs ( variable ) . size ();
  // Return dimension if no out-edges
  if ( num_outedges == 0 ) {
    return network () . size ();
  }

  // If no extra thresholds (no self-edge threshold blowup)
  if ( network () . num_thresholds ( variable ) == num_outedges ) {
    uint64_t out_index = order () [ variable ] ( threshold );
    return network () . outputs ( variable ) [ out_index ];
  }

  // The variable "threshold" correspond to the domain threshold.
  // In case of self-edge threshold blowup this threshold may not
  // be one of the out-edge thresholds. So we find the corresponding
  // out-edge threshold and use it to get the target node.

  // Get self-edge threshold (return dimension if no self-edge)
  auto self_edge_threshold = [&]() {
    for ( uint64_t thresh = 0; thresh < num_outedges; ++ thresh ) {
      uint64_t out_index = order () [ variable ] ( thresh );
      uint64_t target = network () . outputs ( variable ) [ out_index ];
      if ( variable == target ) { // If self-edge
        return thresh; // Return thresh
      }
    }
    // No self-edge (return dimension)
    return network () . size ();
  };

  uint64_t out_edge_threshold = threshold;
  // If threshold is bigger than the self-edge threshold
  // we need to subtract the extra threshold from it
  if ( self_edge_threshold () < threshold ) {
    out_edge_threshold -= 1;
  }

  uint64_t out_index = order () [ variable ] ( out_edge_threshold );
  return network () . outputs ( variable ) [ out_index ];
}

INLINE_IF_HEADER_ONLY std::vector<uint64_t> Parameter::
labelling ( void ) const {
  std::vector<uint64_t> result;
  uint64_t D = network () . size ();

  // pre-allocated vectors (for efficiency)
  std::vector<uint64_t> lower_limits ( D );
  std::vector<uint64_t> upper_limits ( D );
  std::vector<uint64_t> dom ( D );
  std::vector<uint64_t> width ( D );

  std::vector<uint64_t> limits ( D );
  // std::vector<uint64_t> limits = network () . domains ();
  std::vector<uint64_t> jump ( D ); // index offset in each dim
  uint64_t N = 1;
  for ( uint64_t d = 0; d < D; ++ d ) {
    // uint64_t m = network() . outputs ( d ) . size ();
    // Treat the no out edge case as one out edge
    // limits [ d ] = (m ? m : 1) + 1;
    uint64_t m = 1;
    limits [ d ] = m + 1;
    jump[d] =  N;
    N *= limits [ d ];
  }
  // N is now number of domains
  // Domains are implicitly indexed.
  // "jump" is an array telling us how much to change the index
  //   to move +1 in each dimension
  result . resize ( N, 0 );
  for ( int d = 0; d < D; ++ d ) {
    uint64_t n = network() . inputs ( d ) . size ();
    uint64_t numInComb = ( 1LL << n );
    //uint64_t num_debug_domains = 0;
    for ( uint64_t in = 0; in < numInComb; ++ in ) {
      /// What bin does the target point land in for dimension d?
      uint64_t bin = data_ -> logic_ [ d ] . bin ( in );
      /// Which domains have this input combination for dimension d?
      std::fill ( lower_limits.begin(), lower_limits.end(), 0 );
      upper_limits = limits;
      uint64_t sources = network () . inputs ( d ) . size ();
      for ( uint64_t inorder = 0; inorder < sources; ++ inorder ) {
        //std::cout << "Dim " << d << ", in = " << in << " inorder = " << inorder << "\n";
        uint64_t source = network () . inputs ( d ) [ inorder ];
        //std::cout << "source = " << source << "\n";
        // bool activating = network () . interaction ( source, d );
        //std::cout << "high is activating? " << ( activating ? "yes" : "no" ) << "\n";
        // int outorder = network () . order ( source, d );
        //std::cout << "outorder = " << outorder << "\n";
        bool side = in & ( 1LL << inorder );
        //std::cout << "on activating side? " << ( side ? "yes" : "no" ) << "\n";
        // uint64_t thres = data_ -> order_ [ source ] . inverse ( outorder ) + 1;
        uint64_t thres = 1;
        //std::cout << "critical bin = " << thres << "\n";
        // if ( activating ^ side ) {
        if ( not side ) {
          //std::cout << "Case A.\n";
          lower_limits[source] = 0;
          upper_limits[source] = thres;
        } else {
          //std::cout << "Case B.\n";
          lower_limits[source] = thres;
          upper_limits[source] = limits [ source ];
        }
      }
      /// Iterate through two zones:
      ///   Zone 1. domain left of bin
      ///   Zone 2. domain right of bin
      ///   Note. domains matching bin do not
      ///         require anything to be done
      auto apply_mask = [&] ( uint64_t mask ) {
        // Iterate between lower and upper limits applying mask
        uint64_t dom_index = 0;
        dom = lower_limits;
        for ( uint64_t k = 0; k < D; ++ k ) {
          width[k] = upper_limits[k] - lower_limits[k];
          dom_index += jump[k] * lower_limits[k];
          if ( width[k] == 0 ) return;
        }
        while ( 1 ) {
          // apply mask
          result [ dom_index ] |= mask;
          // next domain
          for ( uint64_t k = 0; k < D; ++ k ) {
            ++ dom[k];
            dom_index += jump[k];
            if ( dom[k] < upper_limits[k] ) break;
            dom[k] = lower_limits[k];
            dom_index -= width[k] * jump[k];
          }
          // If we are back to start, return
          bool done = true;
          for ( uint64_t k = 0; k < D; ++ k ) {
            if ( dom[k] != lower_limits[k] ) done = false;
          }
          if ( done ) break;
        }
      };

      uint64_t left = lower_limits [ d ];
      uint64_t right = upper_limits [ d ];

      // Zone 1 (Flows to right)
      if ( bin > left ) {
        lower_limits [ d ] = left;
        // Bug fix for self repressor case
        upper_limits [ d ] = std::min ( right, bin );
        // upper_limits [ d ] = bin;
        apply_mask (1LL << (D+d));
      }
      // Zone 2 (Flows to left)
      if ( bin + 1 < right ) {
        // Bug fix for self repressor case
        lower_limits [ d ] = std::max ( left, bin + 1 );
        // lower_limits [ d ] = bin + 1;
        upper_limits [ d ] = right;
        apply_mask (1LL << d);
      }
    }
  }

  return result;
}

INLINE_IF_HEADER_ONLY Network const Parameter::
network ( void ) const {
  return data_ -> network_;
}

INLINE_IF_HEADER_ONLY std::string Parameter::
stringify ( void ) const {
  std::stringstream ss;
  uint64_t D = data_ -> network_ . size ();
  ss << "[";
  for ( uint64_t d = 0; d < D; ++ d ) {
    if ( d > 0 ) ss << ",";
    ss << "[\"" << network() . name ( d ) << "\","
       << data_ -> logic_[d] << "," << data_ -> order_[d] << "]";
  }
  ss << "]";
  return ss . str ();
}

INLINE_IF_HEADER_ONLY void Parameter::
parse ( std::string const& str ) {
  json p = json::parse(str);
  data_ -> logic_ . clear ();
  data_ -> order_ . clear ();
  for ( auto network_node : p ) {
    data_ -> logic_ . push_back ( LogicParameter () );
    data_ -> order_ . push_back ( OrderParameter () );
    // (*node)[0] has node name, which we ignore.
    data_ -> logic_ . back() . parse ( json::stringify ( network_node[1] )); //TODO inefficient
    data_ -> order_ . back() . parse ( json::stringify ( network_node[2] )); //TODO inefficient
  }
}

INLINE_IF_HEADER_ONLY std::string Parameter::
input_polynomial ( uint64_t in, uint64_t d ) const {
  std::stringstream input_ss;
  std::vector<std::vector<uint64_t>> logic = network () . logic ( d );
  std::string const& node_name = network() . name ( d );
  uint64_t n = network() . inputs ( d ) . size ();
  input_ss << "p" << in << " = ";
  // Corner case: n == 0 (no inputs)
  if ( n == 0 ) {
    input_ss << "B[" << node_name << "]";
    return input_ss . str ();
  }
  uint64_t bit = 1;
  uint64_t k = 0;
  bool first_entry = true;
  for ( auto const& factor : logic ) {
    // if ( factor . size () > 1 ) input_ss << "(";
    bool inner_first = true;
    for ( uint64_t source : factor ) {
      // if ( inner_first ) inner_first = false; else input_ss << " + ";
      std::string source_name = network() . name( source );
      if ( in & bit ) {
        if ( inner_first ) {
          inner_first = false;
          input_ss << (network() . interaction ( source, d ) ? "" : " - ");
        } else {
          input_ss << (network() . interaction ( source, d ) ? " + " : " - ");
        }
        input_ss << "W[" << source_name << "->" << node_name << "]";
        first_entry = false;
      }
      // if ( in & bit ) {
      //   input_ss << "U[" << source_name << "->" << node_name << "]";
      // } else {
      //   input_ss << "L[" << source_name << "->" << node_name << "]";
      // }
      bit <<= 1;
      ++ k;
    }
    // if ( factor . size () > 1 ) input_ss << ")";
    // else if ( k < n ) input_ss << " ";
  }
  input_ss << (first_entry ? "" : " + ") << "b[" << node_name << "]";
  return input_ss . str ();
}

INLINE_IF_HEADER_ONLY std::string Parameter::
output_threshold ( uint64_t j, uint64_t d ) const {
  uint64_t target = network() . outputs ( d ) [ j ];
  std::string node_name = network() . name ( d );
  std::string target_name = network() . name(target);
  std::stringstream output_ss;
  output_ss << "t" << j << " = T[" << node_name << "->" << target_name << "]";
  return output_ss . str ();
}

INLINE_IF_HEADER_ONLY std::string Parameter::
partialorders ( std::string const& type ) const {
  // The default value for type is "", which uses the default type "t"
  std::string thres_type = type.empty() ? "t" : type;
  if ( not ( thres_type == "t" or thres_type == "T" ) ) {
    throw std::runtime_error ( "Invalid threshold type!" );
  }
  // Print parameter partial order
  uint64_t D = data_ -> network_ . size ();
  std::stringstream result_ss;
  for ( uint64_t d = 0; d < D; ++ d ) {
    uint64_t n = network() . inputs ( d ) . size ();
    uint64_t m = 1;
    // uint64_t m = network() . outputs ( d ) . size ();
    uint64_t N = ( 1LL << n );
    // Upper bound threshold for p_i
    std::vector<uint64_t> upper_thres (N);
    // Get order of input combinations and thresholds
    for ( uint64_t i = 0; i < N; ++ i ) {
      uint64_t j = 0;
      while ( j < m && data_ -> logic_[d] ( i * m + j ) ) ++ j;
      upper_thres [i] = j; // Threshold such that p_i < T_j
    }
    // Get partial order as a vector of string
    std::vector<std::string> partial_order;
    for ( uint64_t j = 0; j <= m; ++ j ) {
      for ( uint64_t i = 0; i < N; ++ i ) {
        if ( upper_thres [i] == j ) { // If p_i < T_j
          std::stringstream p_ss;
          p_ss << "p" << i;
          partial_order . push_back ( p_ss . str () );
        }
      }
      if ( j < m ) { // If j == m then p_i > all thresholds
        // Output threshold string
        std::string out_thres;
        if ( network() . outputs ( d ) . size () == 0 ) {
          std::string node_name = network() . name ( d );
          std::stringstream output_ss;
          output_ss << "t" << j << " = T[" << node_name << "->]";
          out_thres = output_ss . str ();
        }
        else {
          // Get original out edge order
          uint64_t j0 = data_ -> order_[d](j);
          // Get output threshold in the format "tj = T[x->y]""
          out_thres = output_threshold ( j0, d );
        }
        // Split string at char '=' to get tj and T[x->y]
        std::stringstream thres_ss (out_thres);
        std::string thres_t_str; // Get tj
        std::getline (thres_ss, thres_t_str, '=');
        std::string thres_T_str; // Get T[x->y]
        std::getline (thres_ss, thres_T_str, '=');
        if ( thres_type == "t" ) {
          // Remove trailing whitespaces
          thres_t_str = std::regex_replace(thres_t_str, std::regex(" +$"), "");
          partial_order . push_back ( thres_t_str ); // Use tj
        } else { // thres_type == "T"
          // Remove leading whitespaces
          thres_T_str = std::regex_replace(thres_T_str, std::regex("^ +"), "");
          partial_order . push_back ( thres_T_str ); // Use T[x->y]
        }
      }
    }
    // Form output string
    std::string node_name = network() . name ( d );
    result_ss << node_name << " : (";
    bool first = true;
    for ( auto item_str : partial_order ) {
      if ( first ) {
        result_ss << item_str;
        first = false;
      } else {
        result_ss << ", " << item_str;
      }
    }
    if ( d < D - 1 ) {
      result_ss << ")\n";  
    } else {
      result_ss << ")";
    }
  }
  return result_ss . str ();
}

INLINE_IF_HEADER_ONLY std::vector<LogicParameter> const & Parameter::
logic ( void ) const {
    return data_ -> logic_;
}

INLINE_IF_HEADER_ONLY std::vector<OrderParameter> const & Parameter::
order ( void ) const {
    return data_ -> order_;
}

INLINE_IF_HEADER_ONLY std::ostream& operator << ( std::ostream& stream, Parameter const& p ) {
  stream << p.stringify();
  return stream;
}
