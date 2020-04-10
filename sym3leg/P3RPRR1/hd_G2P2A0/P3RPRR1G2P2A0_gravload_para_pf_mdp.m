% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G2P2A0
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
%   see P3RPRR1G2P2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G2P2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:17
% EndTime: 2020-03-09 21:25:18
% DurationCPUTime: 0.73s
% Computational Cost: add. (945->146), mult. (669->180), div. (81->4), fcn. (570->63), ass. (0->107)
t2101 = MDP(4) * pkin(1);
t2111 = MDP(2) / 0.2e1 + t2101 / 0.2e1;
t2109 = MDP(3) / 0.2e1;
t2108 = MDP(6) / 0.2e1;
t2107 = MDP(7) / 0.2e1;
t2073 = 0.1e1 / pkin(3);
t2106 = t2073 / 0.2e1;
t2063 = legFrame(3,2);
t2054 = sin(t2063);
t2057 = cos(t2063);
t2009 = t2057 * g(1) - t2054 * g(2);
t2077 = pkin(7) + qJ(3,3);
t2045 = qJ(1,3) + t2077;
t2039 = sin(t2045);
t2042 = cos(t2045);
t1982 = g(3) * t2039 - t2009 * t2042;
t2012 = 0.1e1 / (pkin(1) * sin(t2077) + sin(qJ(3,3)) * pkin(2));
t2100 = t1982 * t2012;
t2064 = legFrame(2,2);
t2055 = sin(t2064);
t2058 = cos(t2064);
t2010 = t2058 * g(1) - t2055 * g(2);
t2078 = pkin(7) + qJ(3,2);
t2046 = qJ(1,2) + t2078;
t2040 = sin(t2046);
t2043 = cos(t2046);
t1983 = g(3) * t2040 - t2010 * t2043;
t2013 = 0.1e1 / (pkin(1) * sin(t2078) + sin(qJ(3,2)) * pkin(2));
t2099 = t1983 * t2013;
t2065 = legFrame(1,2);
t2056 = sin(t2065);
t2059 = cos(t2065);
t2011 = t2059 * g(1) - t2056 * g(2);
t2079 = pkin(7) + qJ(3,1);
t2047 = qJ(1,1) + t2079;
t2041 = sin(t2047);
t2044 = cos(t2047);
t1984 = g(3) * t2041 - t2011 * t2044;
t2014 = 0.1e1 / (pkin(1) * sin(t2079) + sin(qJ(3,1)) * pkin(2));
t2098 = t1984 * t2014;
t1985 = g(3) * t2042 + t2009 * t2039;
t2097 = t1985 * t2012;
t1986 = g(3) * t2043 + t2010 * t2040;
t2096 = t1986 * t2013;
t1987 = g(3) * t2044 + t2011 * t2041;
t2095 = t1987 * t2014;
t2069 = cos(qJ(1,3));
t2080 = qJ(1,3) + pkin(7);
t2094 = (-pkin(3) * t2042 - pkin(2) * cos(t2080) - t2069 * pkin(1)) * t2012;
t2070 = cos(qJ(1,2));
t2081 = qJ(1,2) + pkin(7);
t2093 = (-pkin(3) * t2043 - pkin(2) * cos(t2081) - t2070 * pkin(1)) * t2013;
t2071 = cos(qJ(1,1));
t2082 = qJ(1,1) + pkin(7);
t2092 = (-pkin(3) * t2044 - pkin(2) * cos(t2082) - t2071 * pkin(1)) * t2014;
t2033 = t2063 + t2080;
t2027 = qJ(3,3) + t2033;
t2034 = -t2063 + t2080;
t2028 = qJ(3,3) + t2034;
t2000 = sin(t2027) + sin(t2028);
t2091 = t2000 * t2012;
t2035 = t2064 + t2081;
t2029 = qJ(3,2) + t2035;
t2036 = -t2064 + t2081;
t2030 = qJ(3,2) + t2036;
t2001 = sin(t2029) + sin(t2030);
t2090 = t2001 * t2013;
t2037 = t2065 + t2082;
t2031 = qJ(3,1) + t2037;
t2038 = -t2065 + t2082;
t2032 = qJ(3,1) + t2038;
t2002 = sin(t2031) + sin(t2032);
t2089 = t2002 * t2014;
t2003 = -cos(t2028) + cos(t2027);
t2088 = t2003 * t2012;
t2004 = -cos(t2030) + cos(t2029);
t2087 = t2004 * t2013;
t2005 = -cos(t2032) + cos(t2031);
t2086 = t2005 * t2014;
t2085 = t2012 * t2042;
t2084 = t2013 * t2043;
t2083 = t2014 * t2044;
t2066 = sin(qJ(1,3));
t2076 = g(3) * t2066 - t2009 * t2069;
t2067 = sin(qJ(1,2));
t2075 = g(3) * t2067 - t2010 * t2070;
t2068 = sin(qJ(1,1));
t2074 = g(3) * t2068 - t2011 * t2071;
t2053 = qJ(1,1) - t2065;
t2052 = qJ(1,1) + t2065;
t2051 = qJ(1,2) - t2064;
t2050 = qJ(1,2) + t2064;
t2049 = qJ(1,3) - t2063;
t2048 = qJ(1,3) + t2063;
t2008 = -t2056 * g(1) - t2059 * g(2);
t2007 = -t2055 * g(1) - t2058 * g(2);
t2006 = -t2054 * g(1) - t2057 * g(2);
t1996 = g(3) * t2071 + t2011 * t2068;
t1995 = g(3) * t2070 + t2010 * t2067;
t1994 = g(3) * t2069 + t2009 * t2066;
t1981 = -t2005 * pkin(3) + (cos(t2038) - cos(t2037)) * pkin(2) + (cos(t2053) - cos(t2052)) * pkin(1);
t1980 = -t2004 * pkin(3) + (cos(t2036) - cos(t2035)) * pkin(2) + (cos(t2051) - cos(t2050)) * pkin(1);
t1979 = -t2003 * pkin(3) + (cos(t2034) - cos(t2033)) * pkin(2) + (cos(t2049) - cos(t2048)) * pkin(1);
t1978 = -t2002 * pkin(3) + (-sin(t2037) - sin(t2038)) * pkin(2) + (-sin(t2052) - sin(t2053)) * pkin(1);
t1977 = -t2001 * pkin(3) + (-sin(t2035) - sin(t2036)) * pkin(2) + (-sin(t2050) - sin(t2051)) * pkin(1);
t1976 = -t2000 * pkin(3) + (-sin(t2033) - sin(t2034)) * pkin(2) + (-sin(t2048) - sin(t2049)) * pkin(1);
t1 = [(t2054 * t2006 + t2055 * t2007 + t2056 * t2008) * MDP(4) - g(1) * MDP(8) + (t1994 * t2091 + t1995 * t2090 + t1996 * t2089) * t2109 + (t1982 * t2091 + t1983 * t2090 + t1984 * t2089) * t2108 + (t1985 * t2091 + t1986 * t2090 + t1987 * t2089) * t2107 + ((t1976 * t2100 + t1977 * t2099 + t1978 * t2098) * MDP(6) + (t1976 * t2097 + t1977 * t2096 + t1978 * t2095) * MDP(7)) * t2106 + t2111 * (t2074 * t2089 + t2075 * t2090 + t2076 * t2091); (t2057 * t2006 + t2058 * t2007 + t2059 * t2008) * MDP(4) - g(2) * MDP(8) + (t1994 * t2088 + t1995 * t2087 + t1996 * t2086) * t2109 + (t1982 * t2088 + t1983 * t2087 + t1984 * t2086) * t2108 + (t1985 * t2088 + t1986 * t2087 + t1987 * t2086) * t2107 + ((t1979 * t2100 + t1980 * t2099 + t1981 * t2098) * MDP(6) + (t1979 * t2097 + t1980 * t2096 + t1981 * t2095) * MDP(7)) * t2106 + t2111 * (t2074 * t2086 + t2075 * t2087 + t2076 * t2088); (t1994 * t2085 + t1995 * t2084 + t1996 * t2083) * MDP(3) + (t1982 * t2085 + t1983 * t2084 + t1984 * t2083) * MDP(6) + (t1985 * t2085 + t1986 * t2084 + t1987 * t2083) * MDP(7) - g(3) * MDP(8) + ((t1982 * t2094 + t1983 * t2093 + t1984 * t2092) * MDP(6) + (t1985 * t2094 + t1986 * t2093 + t1987 * t2092) * MDP(7)) * t2073 + (MDP(2) + t2101) * (t2074 * t2083 + t2075 * t2084 + t2076 * t2085);];
taugX  = t1;
