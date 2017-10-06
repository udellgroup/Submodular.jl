#############################################################################
# modularity.jl
# This file handles the basic rules on interactions of mathematical expressions
# to create new expressions concerning their modularity.
#
# For example: negative of a SubModularity expression is SuperModularity.
#
# It's not an all-inclusive file, and set functions whose modularity is not determined
# are classified as NotDetermined.
#############################################################################

import Base: +, -, *, sign

export ExprModularity, SubModularity, SuperModularity, Modularity, ConstModularity, NotDetermined
export +, -, *

# ExprModularity subtypes
abstract type ExprModularity		        end
type SubModularity <: ExprModularity		end
type SuperModularity <: ExprModularity	end
type Modularity <: ExprModularity       end
type ConstModularity <: ExprModularity  end

type NotDetermined <: Vexity
	function NotDetermined()
		warn("Expression's modularity is not determined. Trying to solve problems with uncertain modularity can lead to unexpected behavior.")
    return new()
	end
end

+(m::ExprModularity) = m

+(m::NotDetermined, n::ExprModularity) = m
+(m::ExprModularity, n::NotDetermined) = n
+(m::NotDetermined, n::NotDetermined) = m

+(m::ConstModularity, n::ConstModularity) = m
+(m::ConstModularity, n::SubModularity) = n
+(m::SubModularity, n::ConstModularity) = m
+(m::ConstModularity, n::SuperModularity) = n
+(m::SuperModularity, n::ConstModularity) = m

+(m::Modularity, n::Modularity) = m
+(m::Modularity, n::SubModularity) = n
+(m::SubModularity, n::Modularity) = m
+(m::Modularity, n::SuperModularity) = NotDetermined()
+(m::SuperModularity, n::Modularity) = NotDetermined()

-(m::ExprModularity) = m
-(m::SubModularity) = SuperModularity()
-(m::SuperModularity) = SubModularity()

-(m::ExprModularity, n::ExprModularity) = m + (-n)
