#include <iostream>
#include <vector>
#include <random>
#include <time.h>

/*
You are given a locked container represented as a two-dimensional grid of boolean values (true = locked, false = unlocked).
Your task is to write an algorithm that fully unlocks the box, i.e.,
transforms the entire matrix into all false.

Implement the function:
bool openBox(uint32_t y, uint32_t x);
This function should:
    - Use the SecureBox public API (toggle, isLocked, getState).
    - Strategically toggle cells to reach a state where all elements are false.
    - Return true if the box remains locked, false if successfully unlocked.
You are not allowed to path or modify the SecureBox class.

Evaluation Criteria:
    - Functional correctness
    - Computational efficiency
    - Code quality, structure, and comments
    - Algorithmic insight and clarity
*/

class SecureBox
{
private:
    std::vector<std::vector<bool>> box;

public:

    //================================================================================
    // Constructor: SecureBox
    // Description: Initializes the secure box with a given size and
    //              shuffles its state using a pseudo-random number generator
    //              seeded with current time.
    //================================================================================
    SecureBox(uint32_t y, uint32_t x): ySize(y), xSize(x)
    {
        rng.seed(time(0));
        box.resize(y);
        for (auto& it : box)
            it.resize(x);
        shuffle();
    }

    //================================================================================
    // Method: toggle
    // Description: Toggles the state at position (x, y) and also all cells in the
    //              same row above and the same column to the left of it.
    //================================================================================
    void toggle(uint32_t y, uint32_t x)
    {
        box[y][x] = !box[y][x];
        for (uint32_t i = 0; i < xSize; i++)
            box[y][i] = !box[y][i];
        for (uint32_t i = 0; i < ySize; i++)
            box[i][x] = !box[i][x];
    }

    //================================================================================
    // Method: isLocked
    // Description: Returns true if any cell
    //              in the box is true (locked); false otherwise.
    //================================================================================
    bool isLocked()
    {
        for (uint32_t x = 0; x < xSize; x++)
            for (uint32_t y = 0; y < ySize; y++)
                if (box[y][x])
                    return true;

        return false;
    }

    //================================================================================
    // Method: getState
    // Description: Returns a copy of the current state of the box.
    //================================================================================
    std::vector<std::vector<bool>> getState()
    {
        return box;
    }

private:
    std::mt19937_64 rng;
    uint32_t ySize, xSize;

    //================================================================================
    // Method: shuffle
    // Description: Randomly toggles cells in the box to
    // create an initial locked state.
    //================================================================================
    void shuffle()
    {
        for (uint32_t t = rng() % 1000; t > 0; t--)
            toggle(rng() % ySize, rng() % xSize);
    }
};

//================================================================================
// Function: openBox
// Description: Your task is to implement this function to unlock the SecureBox.
//              Use only the public methods of SecureBox (toggle, getState, isLocked).
//              You must determine the correct sequence of toggle operations to make
//              all values in the box 'false'. The function should return false if
//              the box is successfully unlocked, or true if any cell remains locked.
// Updated Description: Models the toggle operations of the SecureBox as a system of linear equations (mod 2: 1 and 0)
//              and solves it using Gaussian elimination with XOR operations. Each cell's toggle state
//              is represented as a binary variable, and the system incorporates constraints based on
//              the 'toggle' function's behavior, which affects the cell and its entire row and column.
//              The function first applies forward elimination to reduce the system to an upper triangular
//              form, then uses back substitution to determine the sequence of 'toggle' operations required
//              to unlock the box. Finally, it checks whether all cells are unlocked by applying the
//              solution and returns false if the box is unlocked or true if any cell remains locked.
//================================================================================
bool openBox(uint32_t y, uint32_t x) {
    SecureBox box(y, x);
    const auto& state = box.getState();
    
    // System of linear equations: (y * x) x (y * x + 1)
    std::vector<std::vector<int>> equations(y * x, std::vector<int>(y * x + 1, 0));

    for (uint32_t i = 0; i < y; ++i) {
        for (uint32_t j = 0; j < x; ++j) {
            int idx = i * x + j;
            equations[idx][idx] = 1; // always 1 on the diagonal (because each cell toggles itself)
            
            // Button affects all neighbors - in the same row and column
            for (uint32_t k = 0; k < x; ++k) {
                if (k != j) {
                    equations[idx][i * x + k] = 1;
                }
            }
            for (uint32_t k = 0; k < y; ++k) {
                if (k != i) {
                    equations[idx][k * x + j] = 1;
                }
            }
            
            // Recording what we need to get - the cell's state to be open
            equations[idx][y * x] = state[i][j] ? 1 : 0;
        }
    }
    
    // Solving the equation system to find which buttons to press
    for (uint32_t col = 0; col < y * x; ++col) {
        uint32_t row = col;
        
        //Looking for a pivot (key) button
        while (row < y * x && equations[row][col] == 0) ++row;
        if (row == y * x) continue;

        std::swap(equations[row], equations[col]); // Swapping rows to put the pivot in place

        // Deleting the extras to simplify the equations
        for (uint32_t i = col + 1; i < y * x; ++i) {
            if (equations[i][col] == 1) {
                for (uint32_t k = 0; k <= y * x; ++k) {
                    equations[i][k] ^= equations[col][k];
                }
            }
        }
    }

    // Back substitution to derive solution vector (finding which buttons to press to open everything)
    std::vector<bool> solution(y * x, false);
 
    for (int row = y * x - 1; row >= 0; --row) {
        int pivot = -1;
        
        // Locate pivot in row for back substitution
        for (uint32_t col = 0; col < y * x; ++col) {
            if (equations[row][col] == 1) {
                pivot = col;
                break;
            }
        }

        if (pivot == -1) continue; // Skip if row is zero-vector
        solution[pivot] = equations[row][y * x];

        // Eliminate entries above pivot via XOR
        for (int prev = row - 1; prev >= 0; --prev) {
            if (equations[prev][pivot] == 1) {
                equations[prev][y * x] ^= solution[pivot];
            }
        }
    }
    
    // Apply the solution
    for (uint32_t i = 0; i < y; ++i) {
        for (uint32_t j = 0; j < x; ++j) {
            if (solution[i * x + j]) {
                box.toggle(i, j);
            }
        }
    }

    return box.isLocked();
}


int main(int argc, char* argv[])
{
    uint32_t y = std::atol(argv[1]);
    uint32_t x = std::atol(argv[2]);
    bool state = openBox(y, x);

    if (state)
        std::cout << "BOX: LOCKED!" << std::endl;
    else
        std::cout << "BOX: OPENED!" << std::endl;

    return state;
}

