function ranking = RankPop(pop, key)
  for l=1:length(pop);
    fs(l) = pop(l).f;
  end;
  [sfs, ranking] = sort(fs, key);
end
