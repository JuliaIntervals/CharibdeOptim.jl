"""GenerateRandom function can give any random Integer between min and max except values in x """

function GenerateRandom(min, max, x...)

   r = rand(min:max)
   if r in Set(x)
      return generateRandom(min, max, x...)
   else
      return r
end
