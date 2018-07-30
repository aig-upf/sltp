#include <cassert>
#include <iostream>

using namespace std;

class Gripper {
  protected:
    int num_balls_;

  public:
    Gripper(int num_balls) : num_balls_(num_balls) { }
    ~Gripper() { }

    string ball(int i) const {
        assert(i > 0);
        return string("b") + to_string(i);
    }

    void print_move(ostream &os, int from, int to) const {
        // PRE
        os << from << "  - - ";
        for( int i = 0; i < num_balls_; ++i )
            os << " -";

        // arrow
        os << "  =>  ";

        // POST
        os << to << "  - - ";
        for( int i = 0; i < num_balls_; ++i )
            os << " -";
        os << "  ";

        // label
        os << "LABEL move_" << from << "_" << to << endl;
    }

    void print_pick(ostream &os, int r, int b, char gripper) const {
        // PRE: robby
        os << r << "  ";

        // PRE: grippers
        os << (gripper == 'L' ? "0 -" : "- 0");
        os << " ";

        // PRE: balls
        for( int i = 0; i < num_balls_; ++i ) {
            if( i == b )
                os << " " << r;
            else
                os << " -";
        }
                
        // arrow
        os << "  =>  ";

        // POST: robby
        os << r << "  ";

        // POST: grippers
        if( gripper == 'L' )
            os << ball(1 + b) << " -";
        else
            os << "- " << ball(1 + b);
        os << " ";

        // POST: balls
        for( int i = 0; i < num_balls_; ++i ) {
            if( i == b )
                os << " " << gripper;
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL pick_" << r << "_" << ball(1 + b) << "_" << gripper << endl;
    }

    void print_drop(ostream &os, int r, int b, char gripper) const {
        // PRE: robby
        os << r << "  ";

        // PRE: grippers
        if( gripper == 'L' )
            os << ball(1 + b) << " -";
        else
            os << "- " << ball(1 + b);
        os << " ";

        // PRE: balls
        for( int i = 0; i < num_balls_; ++i ) {
            if( i == b )
                os << " " << gripper;
            else
                os << " -";
        }
                
        // arrow
        os << "  =>  ";

        // POST: robby
        os << r << "  ";

        // POST: grippers
        os << (gripper == 'L' ? "0 -" : "- 0");
        os << " ";

        // POST: balls
        for( int i = 0; i < num_balls_; ++i ) {
            if( i == b )
                os << " " << r;
            else
                os << " -";
        }
        os << "  ";

        // label
        os << "LABEL drop_" << r << "_" << ball(1 + b) << "_" << gripper << endl;
    }

    void print_goal(ostream &os) const {
        os << "GOAL  0  0 0 ";
        for( int i = 0; i < num_balls_; ++i )
            os << " 0";
        os << endl;
    }

    void print_psvn(ostream &os) const {
        os << "# written by Blai Bonet" << endl << endl;

        // domains
        os << "DOMAIN B 4 0 1 L R" << endl;
        os << "DOMAIN X " << 1 + num_balls_ << " 0";
        for( int i = 0; i < num_balls_; ++i ) {
            os << " " << ball(1 + i);
        }
        os << endl << endl;

        // state vector
        os << 3 + num_balls_ << " # <robby> <left-gripper> <right-gripper> <pos(ball)>" << endl;
        os << "2 X X";
        for( int i = 0; i < num_balls_; ++i )
            os << " B";
        os << endl << endl;

        // move(from,to)
        os << "# move(from,to)" << endl;
        print_move(os, 0, 1);
        print_move(os, 1, 0);
        os << endl;

        // pick(room,ball,gripper)
        os << "# pick(room,ball,gripper)" << endl;
        for( int r = 0; r < 2; ++r ) {
            for( int b = 0; b < num_balls_; ++b ) {
                print_pick(os, r, b, 'L');
                print_pick(os, r, b, 'R');
            }
        }
        os << endl;

        // drop(room,ball,gripper)
        os << "# drop(room,ball,gripper)" << endl;
        for( int r = 0; r < 2; ++r ) {
            for( int b = 0; b < num_balls_; ++b ) {
                print_drop(os, r, b, 'L');
                print_drop(os, r, b, 'R');
            }
        }
        os << endl;

        // goal
        print_goal(os);
    }
};

int main(int argc, const char **argv) {
    assert(argc > 1);
    int num_balls = atoi(argv[1]);
    Gripper gripper(num_balls);
    gripper.print_psvn(cout);
    return 0;
}

