function seq = SequenceEH(x)
	DurationEH = 7;
    if x == 0
        seq = zeros(1,DurationEH);
    else
        seq = (x/DurationEH:x/DurationEH:x);
    end
end