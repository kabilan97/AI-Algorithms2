from scipy.stats import norm
import scipy.integrate as integrate
import math
import itertools


def payoff(mu, sigma, p, pplus, pmin):
    """
    Payoff function but then in Python
    """
    if pplus + pmin != 0:
        ropt = pplus / (pplus + pmin)
    else:
        ropt = 0
    fun = norm.ppf
    bigphi = 0
    if ropt != 0:
        bigphi, _ = integrate.quad(fun, 0, ropt)
    return p * mu - sigma * (pplus + pmin) * bigphi


"""Important variables"""
N= 50
mu = [10, 100]
sigma = [5, 50]

"Expected payment for agent A"
agent_one = payoff(mu[0], sigma[0], 5, 100, 4.9)
print("Expected payment for agent A")
print(agent_one)
"Expected payment for agent B"
agent_two = payoff(mu[1], sigma[1], 5, 100, 4.9)
print("Expected payment for agent B")
print(agent_two)
"Expected payment for grand coalition"
grand_coal = payoff(48*mu[0]+2*mu[1], math.sqrt(48*(sigma[0]**2)+2*(sigma[1]**2)), 5, 100, 4.9)
print("Expected payment for grand coalition")
print(grand_coal)

"Compute all indexes the agent B can take"
all_pos = list(itertools.combinations(range(50), 2))

def compute_totalmarg(pos_1, pos_2):
    """
    Compute the marginal value for agent A for all different 48 locations given the new positions of agent B
    :param pos_1: index for position of agent B
    :param pos_2: index for position of other agent B
    :return: the sum of the maraginal values of agent A in all different positions
    """
    sum_mu = 0
    sum_sig = 0
    sum_exp = 0
    marginal = 0
    for i in range(50):
        if grand_coal - sum_exp < 0:
            break

        if i == pos_1 or i == pos_2:
            sum_mu += mu[1]
            sum_sig = math.sqrt(sum_sig**2 + sigma[1]**2)
            sum_exp = payoff(sum_mu, sum_sig, 5, 100, 4.9)
            continue

        sum_mu += mu[0]
        sum_sig = math.sqrt(sum_sig ** 2 + sigma[0] ** 2)
        sum_exp_2 = payoff(sum_mu, sum_sig, 5, 100, 4.9)
        marginal += sum_exp_2 - sum_exp
        sum_exp = sum_exp_2

    return marginal


def compute_shapley():
    """
    IMPORTANT NOTE: If you plan to run this function it will take a significant amount of time, about 30-60 min.
    Compute the total marginal value of agent A for all different positions of B
    :return: Compute the Shapley value by dividing the total marginal value by the length of all possible options for B and the 48 positions A can take
    """
    total_value = 0
    for i in all_pos:
        total_value += compute_totalmarg(i[0],i[1])

    return total_value/len(all_pos)/48


agent_one_shap = compute_shapley()
print("Shapley Value for Agent A:")
print(agent_one_shap)

agent_two_shap = (grand_coal - agent_one_shap*48)/2

print("Shapley Value for Agent B:")
print(agent_two_shap)

