"""GenerateRandom function can give any random Integer between min and max except values in x """

function generate_random(min::Int, max::Int, x...)

   r = rand(min:max)
   if r in Set(x)
      return generate_random(min, max, x...)
   else
      return r
   end
end
