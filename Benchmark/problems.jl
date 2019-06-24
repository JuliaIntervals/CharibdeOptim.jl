struct NormalProblem{T}
    func::Function
    domain::T
end

struct DifficultProblem{T}
    func::Function
    domain::T
end

prblm_1 = NormalProblem{IntervalBox}((x,y,z)->x^3+y^2+z, IntervalBox(2..3, 5..6, 4..8))
prblm_2 = NormalProblem{IntervalBox}((x,y,z)->((x^3+y^2)^4+z^2, IntervalBox(2..3, 5..6, 4..8))
prblm_3 = NormalProblem{IntervalBox}((x,y,z)->((x^3+y^2)^4+z^2, IntervalBox(2..3, 5..6, 4..8))
