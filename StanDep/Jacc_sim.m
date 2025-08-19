function similarity = Jacc_sim(A,B)

 similarity = sum(A.*B)/(sum(A+B) - sum(A.*B));
end
