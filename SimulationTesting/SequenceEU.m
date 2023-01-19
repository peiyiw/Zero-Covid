function seq = SequenceEU(x)
	DurationEU = 7;
    if x == 0
        seq = zeros(1,DurationEU);
    else
        seq = (x/DurationEU:x/DurationEU:x);
    end
end