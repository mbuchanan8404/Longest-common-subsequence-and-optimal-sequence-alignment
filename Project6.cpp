/*
* Matthew Buchanan
* CS-340
* Project 6: Longest common subsequence and optimal sequence alignment
*/

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

/* Class to take two strings and compute the longest common sub-sequence and optimal sequence alignment */
class SequenceAlignment
{
public:
	string x, y, xSeq, ySeq, lcs;
	string similarLetters[27] = {"bp", "ck", "cs", "dt", "ey", "gj", "gk", "iy", "kq", "mn", "sz", "vw", "aa", "ae", "ai", "ao", "au", "ee", "ei", "eo", "eu", "ii", "io", "iu", "oo", "ou", "uu"};
	vector<vector<int>> scoreMatrix;
	vector<vector<int>> tracebackMatrix;//0 = left, 1 = diagonal, 2 = above;
	int match, gap, similar, dissimilar, score;

	SequenceAlignment()
	{
	};

	~SequenceAlignment()
	{
	};

	SequenceAlignment(string x, string y, int match, int gap, int similar, int dissimilar)
	{
		this->x = x;
		this->y = y;
		this->match = match;
		this->gap = gap;
		this->similar = similar;
		this->dissimilar = dissimilar;
		createMats();
	};

	/* Create the matrices that will be used for the score and traceback maps */
	void createMats()
	{
		int d1 = x.size() + 1;
		int d2 = y.size() + 1;
		scoreMatrix.resize(d1, vector<int>(d2, 0));
		tracebackMatrix.resize(d1, vector<int>(d2, 0));
	}

	/* Function to determine longest common sub-sequence of two strings */
	void longestCSS()
	{
		for (int i = 1; i < x.length() + 1; i++)
		{
			for (int j = 1; j < y.length() + 1; j++)
			{
				if (x[i - 1] == y[j - 1])
				{
					scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
					tracebackMatrix[i][j] = 1;
				}
				else if (scoreMatrix[i - 1][j] >= scoreMatrix[i][j - 1])
				{
					scoreMatrix[i][j] = scoreMatrix[i - 1][j];
					tracebackMatrix[i][j] = 2;
				}
				else
				{
					scoreMatrix[i][j] = scoreMatrix[i][j - 1];
					tracebackMatrix[i][j] = 0;
				}
			}
		}
		cout << endl << "Longest common subsequence:" << endl;
		int c = 0;
		printCSS(x.size(), y.size(), c);
		cout << endl << "Length of lcs: " << c << endl;
	}

	/* Function to output the longest common sub-sequence */
	void printCSS(int i, int j, int &counter)
	{
		if (i == 0 || j == 0)
		{
			return;
		}
		if (tracebackMatrix[i][j] == 1)
		{
			printCSS(i - 1, j - 1, counter);
			cout << x[i - 1];
			counter++;
		}
		else if (tracebackMatrix[i][j] == 2)
			printCSS(i - 1, j, counter);
		else printCSS(i, j - 1, counter);
	}

	/* Function to find the optimal sequence alignment for two strings */
	void optimalSA()
	{
		/* Initialize matrix */
		for (int i = 1; i < x.length() + 1; i++)
		{
			scoreMatrix[i][0] = scoreMatrix[i - 1][0] + gap;
			tracebackMatrix[i][0] = 0;// redundant
		}
		for (int i = 1; i < y.length() + 1; i++)
		{
			scoreMatrix[0][i] = scoreMatrix[0][i - 1] + gap;
			tracebackMatrix[0][i] = 2;
		}
		int diagonal, left, above, max, tracebackDir;
		diagonal = left = above = max = tracebackDir = 0;
		for (int i = 1; i < x.length() + 1; i++)
		{
			for (int j = 1; j < y.length() + 1; j++)
			{
				if (x[i - 1] == y[j - 1])//from diagonal
					diagonal = match + scoreMatrix[i - 1][j - 1];//match
				else
				{
					bool found = false;
					for (int k = 0; k < 27 && !found; k++)
					{
						if ((x[i - 1] == similarLetters[k][0] && y[j - 1] == similarLetters[k][1]) || (y[j - 1] == similarLetters[k][0] && x[i - 1] == similarLetters[k][1]))
						{
							diagonal = similar + scoreMatrix[i - 1][j - 1];//similar
							found = true;
						}
					}
					if (found == false)
						diagonal = dissimilar + scoreMatrix[i - 1][j - 1];//dissimilar
				}
				max = diagonal;
				tracebackDir = 1;
				left = scoreMatrix[i - 1][j] + gap;//from left
				if (left > max)
				{
					max = left;
					tracebackDir = 0;
				}
				above = scoreMatrix[i][j - 1] + gap;//from above
				if (above > max)
				{
					max = above;
					tracebackDir = 2;
				}
				scoreMatrix[i][j] = max;
				tracebackMatrix[i][j] = tracebackDir;
			}
		}
	}

	/* Function to perform traceback and output optimal alignment, score, and maps */
	void printOSA()
	{
		/* Perform traceback */
		xSeq = ySeq = "";
		int i = x.size(), j = y.size();
		score = scoreMatrix[i][j];
		int * tracer = &tracebackMatrix[i][j];
		while (tracer != &tracebackMatrix[0][0])
		{
			if (*tracer == 0)
			{
				xSeq = x[i - 1] + xSeq;
				ySeq = '-' + ySeq;
				i--;
			}
			else if (*tracer == 2)
			{
				xSeq = '-' + xSeq;
				ySeq = y[j - 1] + ySeq;
				j--;
			}
			else
			{
				xSeq = x[i - 1] + xSeq;
				ySeq = y[j - 1] + ySeq;
				i--;
				j--;
			}
			tracer = &tracebackMatrix[i][j];
		}
		cout << "Score Matrix" << endl;
		cout << "     " << setw(3) << '-' << "  ";
		for (int i = 0; i < x.size(); i++)
			cout << setw(3) << x[i] << "  ";
		cout << endl << "  -";
		for (int i = 0; i < y.length() + 1; i++)
		{
			if (i > 0)
				cout << setw(3) << y[i - 1] << "  ";
			else cout << "  ";
			for (int j = 0; j < x.length() + 1; j++)
				cout << setw(3) << scoreMatrix[j][i] << "  ";
			cout << endl;
		}
		cout << endl << endl;
		cout << "Traceback Matrix 0 = left, 1 = diagonal, 2 = above" << endl << "     " << setw(3) << '-' << "  ";
		for (int i = 0; i < x.size(); i++)
			cout << setw(3) << x[i] << "  ";
		cout << endl << "  -";
		for (int i = 0; i < y.length() + 1; i++)
		{
			if (i > 0)
				cout << setw(3) << y[i - 1] << "  ";
			else cout << "  ";
			for (int j = 0; j < x.length() + 1; j++)
				cout << setw(3) << tracebackMatrix[j][i] << "  ";
			cout << endl;
		}
		cout << endl << endl;
		cout << "Optimal alignment:" << endl << xSeq << endl;
		cout << ySeq << endl;
		cout << "Alignment score: " << score << endl << endl;
	}
};

/* Function prototypes */
int initialPrompt();
string retrieveString();
void retrieveParameters(int &m, int &g, int &s, int &d);

/* Begin main */
int main()
{
	string firstString, secondString;
	int match, gap, sim, dsim;
	int selection = initialPrompt();
	cout << "For the first string:" << endl;
	firstString = retrieveString();
	cout << "For the second string:" << endl;
	secondString = retrieveString();
	if (selection == 1)
	{
		SequenceAlignment c(firstString, secondString, 0, 0, 0, 0);
		c.longestCSS();
	}
	else
	{
		retrieveParameters(match, gap, sim, dsim);
		SequenceAlignment s(firstString, secondString, match, gap, sim, dsim);
		s.optimalSA();
		s.printOSA();
	}
	system("PAUSE");
	return 0;
}
/* End main */

/* function to retrieve user's choice of LCS or OSA */
int initialPrompt()
{
	int choice = 0;
	while (choice != 1 && choice != 2 && isdigit(choice))
	{
		cout << "Enter 1 to find the longest common subsequence, or 2 to find the optimal sequence alignment:" << endl;
		cin >> choice;
		if (choice != 1 && choice != 2 && isdigit(choice));
			cout << "Invalid number entered. Please try again" << endl;
		cout << endl;
	}
	return choice;
}

/* Retrieve and test string from user */
string retrieveString()
{
	bool correctInput = false;
	string s;
	while (!correctInput)
	{
		correctInput = true;
		s = "";
		cout << "Please enter a string of english language letters:" << endl;
		cin >> s;
		if (s == "")
		{
			correctInput = false;
			cout << "Empty strings are not valid. Please try again:" << endl;//not needed?
		}
		else
		{
			bool nonLetter = false;
			for (int i = 0; i < s.length(); i++)
			{
				if (!isalpha(s[i]))
					nonLetter = true;
			}
			if (nonLetter)
			{
				cout << "Only english language letters are valid. Please try again:" << endl;
				correctInput = false;
			}
		}
	}
	for (int i = 0; i < s.length(); i++)
		s[i] = tolower(s[i]);
	return s;
}

/* Function to retrieve the parameters used for the OSA algorithm */
void retrieveParameters(int &m, int &g, int &s, int &d)
{
	int matchScore, gapPenalty, similarScore, dissimilarPenalty;
	cout << "Please enter a match score of at least 2:" << endl;
	cin >> matchScore;
	while (matchScore < 2)
	{
		cout << "Match score must be greater than 2. Please try again:" << endl;
		cin >> matchScore;
	}
	m = matchScore;
	cout << "Please enter a score for similar letters:" << endl;
	cin >> similarScore;
	while (similarScore >= matchScore || similarScore <= 0)
	{
		cout << "The similar score must be less than the match score and greater than 0. Please try again:" << endl;
		cin >> similarScore;
	}
	s = similarScore;
	cout << "Please enter a penalty for dissimilar letters as a negative integer:" << endl;
	cin >> dissimilarPenalty;
	while (dissimilarPenalty >= 0)
	{
		cout << "The dissimilar penalty must be less than 0. Please try again:" << endl;
		cin >> dissimilarPenalty;
	}
	d = dissimilarPenalty;
	cout << "Finally, please enter a gap penalty as a negative integer:" << endl;
	cin >> gapPenalty;
	while (gapPenalty >= 0)
	{
		cout << "The gap penalty must be less than 0. Please try again:" << endl;
		cin >> gapPenalty;
	}
	cout << endl;
	g = gapPenalty;
}