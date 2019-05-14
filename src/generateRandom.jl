function generateRandom(min, max, x...)
   r = rand(min:max)
   if r in Set(x)
      return generateRandom(min, max, x...)
   else
      return r
end
