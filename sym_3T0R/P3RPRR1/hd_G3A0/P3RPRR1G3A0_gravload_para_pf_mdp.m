% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:15
% EndTime: 2020-03-09 21:27:16
% DurationCPUTime: 1.05s
% Computational Cost: add. (945->146), mult. (669->180), div. (81->4), fcn. (570->63), ass. (0->107)
t2137 = MDP(4) * pkin(1);
t2150 = MDP(2) / 0.2e1 + t2137 / 0.2e1;
t2145 = MDP(3) / 0.2e1;
t2144 = MDP(6) / 0.2e1;
t2143 = MDP(7) / 0.2e1;
t2105 = 0.1e1 / pkin(3);
t2142 = t2105 / 0.2e1;
t2095 = legFrame(3,2);
t2087 = sin(t2095);
t2090 = cos(t2095);
t2042 = t2090 * g(1) - t2087 * g(2);
t2098 = sin(qJ(1,3));
t2101 = cos(qJ(1,3));
t2024 = g(3) * t2101 + t2042 * t2098;
t2108 = pkin(7) + qJ(3,3);
t2078 = qJ(1,3) + t2108;
t2072 = sin(t2078);
t2075 = cos(t2078);
t2015 = g(3) * t2075 + t2042 * t2072;
t2045 = 0.1e1 / (pkin(1) * sin(t2108) + sin(qJ(3,3)) * pkin(2));
t2135 = t2015 * t2045;
t2016 = -g(3) * t2072 + t2042 * t2075;
t2134 = t2016 * t2045;
t2096 = legFrame(2,2);
t2088 = sin(t2096);
t2091 = cos(t2096);
t2043 = t2091 * g(1) - t2088 * g(2);
t2109 = pkin(7) + qJ(3,2);
t2079 = qJ(1,2) + t2109;
t2073 = sin(t2079);
t2076 = cos(t2079);
t2017 = g(3) * t2076 + t2043 * t2073;
t2046 = 0.1e1 / (pkin(1) * sin(t2109) + sin(qJ(3,2)) * pkin(2));
t2133 = t2017 * t2046;
t2018 = -g(3) * t2073 + t2043 * t2076;
t2132 = t2018 * t2046;
t2097 = legFrame(1,2);
t2089 = sin(t2097);
t2092 = cos(t2097);
t2044 = t2092 * g(1) - t2089 * g(2);
t2110 = pkin(7) + qJ(3,1);
t2080 = qJ(1,1) + t2110;
t2074 = sin(t2080);
t2077 = cos(t2080);
t2019 = g(3) * t2077 + t2044 * t2074;
t2047 = 0.1e1 / (pkin(1) * sin(t2110) + sin(qJ(3,1)) * pkin(2));
t2131 = t2019 * t2047;
t2020 = -g(3) * t2074 + t2044 * t2077;
t2130 = t2020 * t2047;
t2114 = qJ(1,3) + pkin(7);
t2129 = (pkin(3) * t2072 + pkin(2) * sin(t2114) + t2098 * pkin(1)) * t2045;
t2099 = sin(qJ(1,2));
t2115 = qJ(1,2) + pkin(7);
t2128 = (pkin(3) * t2073 + pkin(2) * sin(t2115) + t2099 * pkin(1)) * t2046;
t2100 = sin(qJ(1,1));
t2116 = qJ(1,1) + pkin(7);
t2127 = (pkin(3) * t2074 + pkin(2) * sin(t2116) + t2100 * pkin(1)) * t2047;
t2066 = t2095 + t2114;
t2060 = qJ(3,3) + t2066;
t2067 = -t2095 + t2114;
t2061 = qJ(3,3) + t2067;
t2113 = sin(t2060) - sin(t2061);
t2126 = t2113 * t2045;
t2068 = t2096 + t2115;
t2062 = qJ(3,2) + t2068;
t2069 = -t2096 + t2115;
t2063 = qJ(3,2) + t2069;
t2112 = sin(t2062) - sin(t2063);
t2125 = t2112 * t2046;
t2070 = t2097 + t2116;
t2064 = qJ(3,1) + t2070;
t2071 = -t2097 + t2116;
t2065 = qJ(3,1) + t2071;
t2111 = sin(t2064) - sin(t2065);
t2124 = t2111 * t2047;
t2036 = cos(t2061) + cos(t2060);
t2123 = t2036 * t2045;
t2037 = cos(t2063) + cos(t2062);
t2122 = t2037 * t2046;
t2038 = cos(t2065) + cos(t2064);
t2121 = t2038 * t2047;
t2119 = t2045 * t2072;
t2118 = t2046 * t2073;
t2117 = t2047 * t2074;
t2102 = cos(qJ(1,2));
t2107 = g(3) * t2102 + t2043 * t2099;
t2103 = cos(qJ(1,1));
t2106 = g(3) * t2103 + t2044 * t2100;
t2086 = qJ(1,1) - t2097;
t2085 = qJ(1,1) + t2097;
t2084 = qJ(1,2) - t2096;
t2083 = qJ(1,2) + t2096;
t2082 = qJ(1,3) - t2095;
t2081 = qJ(1,3) + t2095;
t2041 = -t2089 * g(1) - t2092 * g(2);
t2040 = -t2088 * g(1) - t2091 * g(2);
t2039 = -t2087 * g(1) - t2090 * g(2);
t2029 = -g(3) * t2100 + t2044 * t2103;
t2027 = -g(3) * t2099 + t2043 * t2102;
t2025 = -g(3) * t2098 + t2042 * t2101;
t2014 = -t2038 * pkin(3) + (-cos(t2071) - cos(t2070)) * pkin(2) + (-cos(t2086) - cos(t2085)) * pkin(1);
t2013 = -t2037 * pkin(3) + (-cos(t2069) - cos(t2068)) * pkin(2) + (-cos(t2084) - cos(t2083)) * pkin(1);
t2012 = -t2036 * pkin(3) + (-cos(t2067) - cos(t2066)) * pkin(2) + (-cos(t2082) - cos(t2081)) * pkin(1);
t2011 = t2111 * pkin(3) + (sin(t2070) - sin(t2071)) * pkin(2) + (sin(t2085) - sin(t2086)) * pkin(1);
t2010 = t2112 * pkin(3) + (sin(t2068) - sin(t2069)) * pkin(2) + (sin(t2083) - sin(t2084)) * pkin(1);
t2009 = t2113 * pkin(3) + (sin(t2066) - sin(t2067)) * pkin(2) + (sin(t2081) - sin(t2082)) * pkin(1);
t1 = [(t2087 * t2039 + t2088 * t2040 + t2089 * t2041) * MDP(4) - g(1) * MDP(8) + (t2025 * t2123 + t2027 * t2122 + t2029 * t2121) * t2145 + (t2015 * t2123 + t2017 * t2122 + t2019 * t2121) * t2144 + (t2016 * t2123 + t2018 * t2122 + t2020 * t2121) * t2143 + ((t2012 * t2135 + t2013 * t2133 + t2014 * t2131) * MDP(6) + (t2012 * t2134 + t2013 * t2132 + t2014 * t2130) * MDP(7)) * t2142 + t2150 * (t2024 * t2123 + t2106 * t2121 + t2107 * t2122); (t2090 * t2039 + t2091 * t2040 + t2092 * t2041) * MDP(4) - g(2) * MDP(8) + (-t2025 * t2126 - t2027 * t2125 - t2029 * t2124) * t2145 + (-t2015 * t2126 - t2017 * t2125 - t2019 * t2124) * t2144 + (-t2016 * t2126 - t2018 * t2125 - t2020 * t2124) * t2143 + ((t2009 * t2135 + t2010 * t2133 + t2011 * t2131) * MDP(6) + (t2009 * t2134 + t2010 * t2132 + t2011 * t2130) * MDP(7)) * t2142 + t2150 * (-t2024 * t2126 - t2106 * t2124 - t2107 * t2125); (-t2025 * t2119 - t2027 * t2118 - t2029 * t2117) * MDP(3) + (-t2015 * t2119 - t2017 * t2118 - t2019 * t2117) * MDP(6) + (-t2016 * t2119 - t2018 * t2118 - t2020 * t2117) * MDP(7) - g(3) * MDP(8) + ((t2015 * t2129 + t2017 * t2128 + t2019 * t2127) * MDP(6) + (t2016 * t2129 + t2018 * t2128 + t2020 * t2127) * MDP(7)) * t2105 + (MDP(2) + t2137) * (-t2024 * t2119 - t2106 * t2117 - t2107 * t2118);];
taugX  = t1;
