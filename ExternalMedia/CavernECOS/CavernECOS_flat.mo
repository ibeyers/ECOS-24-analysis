within ExternalMedia.CavernECOS;

model CavernECOS_flat
  //####################  Imports  #######################
  package H2CoolProp "CoolProp model of H2"
    extends ExternalMedia.Media.CoolPropMedium(mediumName = "Hydrogen", substanceNames = {"H2"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, SpecificEnthalpy(start = 1400));
  end H2CoolProp;  
  import Modelica.SIunits.SpecificEnergy;
  import Modelica.SIunits.Volume;
  import Modelica.SIunits.Energy;
  import Modelica.SIunits.Area;
  import Modelica.SIunits.Length;
  import Modelica.SIunits.Diameter;
  import Modelica.SIunits.Radius;
  import Modelica.SIunits.Density;
  import Modelica.SIunits.MassFlowRate;
  import Modelica.SIunits.Temperature;
  import Modelica.SIunits.Pressure;
  import Modelica.SIunits.Mass;
  import Modelica.SIunits.HeatFlowRate;
  import Modelica.SIunits.CoefficientOfHeatTransfer;
  import Modelica.Constants.g_n;
  import Modelica.Constants.R;
  import Modelica.SIunits.Conversions.from_bar;
  import Modelica.SIunits.Conversions.from_degC;
  import Modelica.SIunits.ThermalDiffusivity;
  import Modelica.SIunits.ThermalConductivity;
  //####################  Parameters  ####################
  parameter SpecificEnergy LHV_hydrogen = 120*10^6 "Lower heating value of hydrogen";
  parameter Real molar_mass_H2 = 0.002016 "Molar Mass needed for mass initialisation";
  parameter Volume V_cavern = 500000 "Volume of cavern";
  parameter Area A_eff = 249000 "Effective surface area of cavern";
  parameter Pressure p_max_op = from_bar(185) "Chosen maximum operating pressure level";
  parameter Pressure p_min_op = from_bar(60) "Chosen minimum operating pressure level";
  // Initial values
  parameter Pressure p_cavern_initial = from_bar(160) "Initial cavern pressure";
  parameter Temperature T_cavern_initial = 320.63 "initial cavern temperature";
  parameter Temperature T_cavern_infinite = 320.63 "infinite temperature of surrounding salt-rock";
  parameter Diameter d_cavern = 46.06588659617807 "cavern diameter";
  parameter Radius r_cavern = d_cavern/2 "cavern radius";
  // heat transfer to surrounding salt rock
  parameter Integer n_steps_saltrock = 71 "Number of increments for the discretization of the cavern wall";
  parameter Length delta_r_saltrock = 0.1 "Radial length of each increment";
  //parameter Temperature T_salt_init[n_steps_saltrock] = ones(n_steps_saltrock)* T_cavern_infinite;
  parameter Temperature T_salt_init[n_steps_saltrock] = {320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63, 320.63};
  parameter Temperature T_wall_init = 320.63 "initial wall temperature";
  parameter ThermalDiffusivity k_saltrock_diff = 3*10^(0 - 6) "thermal diffusivity of salt rock";
  parameter ThermalConductivity k_saltrock_cond = 6 "Thermal Conductivity of salt rock";
  
  //####################  Interfaces/Inputs ######################
  input MassFlowRate m_dot_H2_in(start = 0) "Mass flow rate of hydrogen injected into the cavern";
  input MassFlowRate m_dot_H2_out(start = 1.5) "Mass flow rate of hydrogen withdrawn rom the cavern";
  parameter Temperature T_in = from_degC(40) "Temperature of the hydrogen gas that is injected into the cavern";
  
  //####################  Initial value problem ######################
  parameter Mass m_H2_start = p_cavern_initial*V_cavern*molar_mass_H2/(R*T_cavern_initial) "start value for hydrogen mass in cavern";
  Mass m_H2_initial "Initial hydrogen mass in cavern";
  parameter Density rho_start = m_H2_start/V_cavern "start value for density of hydrogen gas in cavern";
  
  //####################  Variables ######################
  Density rho(start = rho_start);
  Pressure p_cavern(start = p_cavern_initial) "Pressure inside of the cavern";
  Energy E_cavern(displayUnit = "GWh") "Energy inside of the cavern in the form of hydrogen";
  Mass m_H2_cavern(start = m_H2_start) "Mass of hydrogen gas inside of the cavern";
  Temperature T_H2_cavern "Temperature of the hydrogen inside of the cavern";
  HeatFlowRate Q_dot_rock(displayUnit = "MW") "Heat flux exchanged with the surrooundig salt rock";
  H2CoolProp.ThermodynamicState HydrogenState "Thermodynamic state of H2 in cavern";
  H2CoolProp.ThermodynamicState HydrogenStateInput "Thermodynamic state of H2 coming into cavern";
  H2CoolProp.ThermodynamicState HydrogenStateStart "Initial thermodynamic state of H2 in cavern";
  //heat transfer to surrounding saltrock
  Temperature T_salt[n_steps_saltrock](start = T_salt_init) "salt dome temperatures";
  Temperature T_wall(start = T_wall_init) "Uniform temperature of the cavern wall";
  Real Bi = U*delta_r_saltrock/k_saltrock_cond "Biot number";
  CoefficientOfHeatTransfer U "Overall heat transfer coefficient between the hydrogen gas and the cavern wall";
initial equation
  T_H2_cavern = T_cavern_initial;
  m_H2_cavern = m_H2_initial;
//####################  Equation ######################
equation
//cavern states
  HydrogenStateInput = H2CoolProp.setState_pT(p_cavern, T_in);
  HydrogenState = H2CoolProp.setState_dT(rho, T_H2_cavern);
  HydrogenStateStart = H2CoolProp.setState_pT(p_cavern_initial, T_cavern_initial);
  m_H2_initial = HydrogenStateStart.d*V_cavern;
  p_cavern = HydrogenState.p;
  E_cavern = m_H2_cavern*LHV_hydrogen;
//energy in cavern (including cushion gas)
//energy balance
  m_H2_cavern*HydrogenState.cv*der(T_H2_cavern) = Q_dot_rock + m_dot_H2_in*HydrogenStateInput.h - m_dot_H2_out*HydrogenState.h - (HydrogenState.h - p_cavern*1/HydrogenState.d)*der(m_H2_cavern);
//heat transfer cavern-surrounding saltrock
  Q_dot_rock = U*A_eff*(T_wall - T_H2_cavern);
  U = 0.1*((HydrogenState.beta*g_n*(HydrogenState.d)^2*abs(T_wall - T_H2_cavern)*HydrogenState.cp*(HydrogenState.lambda)^2)/HydrogenState.eta)^(1/3);
  T_H2_cavern*Bi = T_salt[1]*((Bi/2) + 1) + T_salt[2]*((Bi/2) - 1);
  for i in 2:n_steps_saltrock - 1 loop
    der(T_salt[i]) = k_saltrock_diff/(delta_r_saltrock^2)*(T_salt[i + 1]*((delta_r_saltrock/(2*(r_cavern + delta_r_saltrock*(i - 1)))) + 1) + T_salt[i - 1]*(1 - (delta_r_saltrock/(2*(r_cavern + delta_r_saltrock*(i - 1))))) - 2*T_salt[i]);
//i-1 because there is no entry at [0]
  end for;
  T_salt[end] = T_cavern_infinite;
// calculation of wall temperature
  T_wall = (T_salt[2] + T_salt[1])/2;
//mass balance
  der(rho)*V_cavern = m_dot_H2_in - m_dot_H2_out;
  m_H2_cavern = rho*V_cavern;
  
//####################  Component Management ######################
//    if p_cavern > p_max_op then
//  terminate("Error Cavern: p_cavern > p_cavern_max!");
// end if;
// if p_cavern < p_min_op then
//  terminate("Error Cavern: p_cavern < p_cavern_min!");
//  end if;
end CavernECOS_flat;
