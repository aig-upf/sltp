#include <cassert>
#include <iostream>

using namespace std;

class Blocks {
  protected:
    int num_blocks_;

  public:
    Blocks(int num_blocks) : num_blocks_(num_blocks) { }
    ~Blocks() { }

    string block(int i) const {
        assert(i >= 0);
        return i == 0 ? string("0") : string("") + char(int('A' + i - 1));
    }

    void print_stack(ostream &os, int x, int y) const {
        // PRE: holding
        os << block(1 + x) << "  ";

        // PRE: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( (i == x) || (i == y) )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // PRE: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }

        // arrow + POST: holding
        os << "  =>  0  ";

        // POST: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else if( i == y )
                os << " 0";
            else
                os << " -";
        }
        os << "  ";

        // POST: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " " << block(1 + y);
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL stack_" << block(1 + x) << "_" << block(1 + y) << endl;
    }

    void print_unstack(ostream &os, int x, int y) const {
        // PRE: holding
        os << "0  ";

        // PRE: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else if( i == y )
                os << " 0";
            else
                os << " -";
        }
        os << "  ";

        // PRE: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " " << block(1 + y);
            else
                os << " -";
        }

        // arrow + POST: holding
        os << "  =>  " << block(1 + x) << "  ";

        // POST: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( (i == x) || (i == y) )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // POST: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL unstack_" << block(1 + x) << "_" << block(1 + y) << endl;
    }

    void print_pick(ostream &os, int x) const {
        // PRE: holding
        os << "0  ";

        // PRE: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // PRE: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }

        // arrow + POST: holding
        os << "  =>  " << block(1 + x) << "  ";

        // POST: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // POST: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL pick_" << block(1 + x) << endl;
    }

    void print_putdown(ostream &os, int x) const {
        // PRE: holding
        os << block(1 + x) << "  ";

        // PRE: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // PRE: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }

        // arrow + POST: holding
        os << "  =>  0  ";

        // POST: clear
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 1";
            else
                os << " -";
        }
        os << "  ";

        // POST: on
        for( int i = 0; i < num_blocks_; ++i ) {
            if( i == x )
                os << " 0";
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL putdown_" << block(1 + x) << endl;
    }

    void print_goal(ostream &os) const {
        os << "GOAL  0  ";
        for( int i = 0; i < num_blocks_; ++i )
            os << " 1";
        os << "  ";
        for( int i = 0; i < num_blocks_; ++i )
            os << " 0";
        os << endl;
    }

    void print_psvn(ostream &os) const {
        os << "# written by Blai Bonet" << endl << endl;

        // domain
        os << "DOMAIN X " << 1 + num_blocks_;
        for( int i = 0; i <= num_blocks_; ++i ) {
            os << " " << block(i);
        }
        os << endl << endl;

        // state vector
        os << 1 + 2 * num_blocks_ << " # <hold> <clear(x)> <on(x)>" << endl;
        os << "X  ";
        for( int i = 0; i < num_blocks_; ++i )
            os << " 2";
        os << "  ";
        for( int i = 0; i < num_blocks_; ++i )
            os << " X";
        os << endl << endl;

        // stack(x,y)
        os << "# Stack(x,y)" << endl;
        for( int x = 0; x < num_blocks_; ++x ) {
            for( int y = 0; y < num_blocks_; ++y ) {
                if( x != y )
                    print_stack(os, x, y);
            }
            os << endl;
        }

        // unstack(x,y)
        os << "# Unstack(x,y)" << endl;
        for( int x = 0; x < num_blocks_; ++x ) {
            for( int y = 0; y < num_blocks_; ++y ) {
                if( x != y )
                    print_unstack(os, x, y);
            }
            os << endl;
        }

        // pick(x)
        os << "# Pick(x)" << endl;
        for( int x = 0; x < num_blocks_; ++x )
            print_pick(os, x);
        os << endl;

        // putdown(x)
        os << "# Putdown(x)" << endl;
        for( int x = 0; x < num_blocks_; ++x )
            print_putdown(os, x);
        os << endl;

        // goal
        print_goal(os);
    }
};

int main(int argc, const char **argv) {
    assert(argc > 1);
    int num_blocks = atoi(argv[1]);
    Blocks blocks(num_blocks);
    blocks.print_psvn(cout);
    return 0;
}

