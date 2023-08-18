% As you can see, this is a hard coded version because I cant
% think of an efficient algorithm to generate a codebook.
% I might create a brute force algorithm to create a codebook for any
% number of transmit antennas

% Each cell within a cell represents the codebook for the given antenna
% configuration

% USE EXAMPLE: mapping(3){1}{2}; Nt=3, first codebook, second pair
function mapping = codebook_gen()
    key=[3,4,8,2]
    % Value for 3
    value{1} = {
		    {[1,2]},
		    {[2,3]}
	    };
    % Value for 4
    value{2}={
		    {[1,2], [3,4]},
		    {[2,3], [1,4]}
	    };
    % Value for 8
    value{3}={
		    {[1,2],[3,4],[5,6],[7,8]},
		    {[2,3],[4,5],[6,7],[1,8]},
		    {[1,5],[2,6],[3,7],[4,8]},
            {[1,3],[2,4],[5,7],[6,8]}
	    };
    % Value for Nt=2; Existing only for debugging purposes
    value{4} = {
		    {[1,2]}
	    };
    mapping=containers.Map(key,value);
end