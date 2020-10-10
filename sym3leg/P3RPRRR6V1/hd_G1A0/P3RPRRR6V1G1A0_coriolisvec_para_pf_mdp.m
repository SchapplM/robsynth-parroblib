% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR6V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:34
% EndTime: 2020-08-06 18:32:39
% DurationCPUTime: 4.84s
% Computational Cost: add. (22005->343), mult. (19983->599), div. (1515->20), fcn. (9759->68), ass. (0->270)
t2107 = xDP(2);
t2108 = xDP(1);
t2109 = (pkin(6) + pkin(5));
t2011 = -2 * pkin(2) * t2107 - 2 * t2108 * t2109;
t2018 = -pkin(2) * t2108 + t2107 * t2109;
t2071 = qJ(1,3) + legFrame(3,3);
t2054 = pkin(7) + t2071;
t2034 = sin(t2054);
t2057 = sin(t2071);
t2060 = cos(t2071);
t2041 = qJ(3,3) + t2054;
t2042 = -qJ(3,3) + t2054;
t2205 = cos(t2041) + cos(t2042);
t2208 = sin(t2041) + sin(t2042);
t2037 = cos(t2054);
t2272 = 0.2e1 * t2037;
t2275 = -2 * pkin(1);
t1971 = t2018 * t2272 + t2011 * t2034 + (t2057 * t2107 + t2060 * t2108) * t2275 + (-t2107 * t2208 - t2108 * t2205) * pkin(3);
t2081 = cos(pkin(7)) * pkin(1);
t2103 = cos(qJ(3,3));
t2266 = t2103 * pkin(3) + pkin(2);
t2019 = t2081 + t2266;
t2012 = 0.1e1 / t2019;
t2070 = t2081 + pkin(2);
t2100 = sin(qJ(3,3));
t2080 = sin(pkin(7)) * pkin(1);
t2069 = t2080 + pkin(5);
t2114 = 0.1e1 / pkin(3);
t2110 = 0.2e1 * qJ(3,3);
t2085 = sin(t2110);
t2074 = pkin(3) * t2085;
t2287 = 2 * pkin(2);
t2004 = t2100 * t2287 + t2074 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1);
t1992 = 0.1e1 / t2004;
t2244 = t1971 * t1992;
t2178 = t2114 * t2244;
t2146 = t2069 * t2178;
t2001 = -t2034 * t2108 + t2037 * t2107;
t2235 = t2001 * t2012;
t2170 = t2100 * t2235;
t2213 = t2100 * t2103;
t2217 = t2070 * t2103;
t2013 = 0.1e1 / t2019 ^ 2;
t2234 = t2001 * t2013;
t2091 = t2103 ^ 2;
t2282 = MDP(6) * (0.2e1 * t2091 - 0.1e1);
t2290 = 0.2e1 * ((MDP(5) * t2213 + t2282) * t2234 - (MDP(10) * (t2103 * t2146 / 0.2e1 + t2070 * t2170) + MDP(11) * (t2217 * t2235 - t2100 * t2146 / 0.2e1)) * t2012) * t1971;
t2072 = qJ(1,2) + legFrame(2,3);
t2055 = pkin(7) + t2072;
t2035 = sin(t2055);
t2058 = sin(t2072);
t2061 = cos(t2072);
t2045 = qJ(3,2) + t2055;
t2046 = -qJ(3,2) + t2055;
t2204 = cos(t2045) + cos(t2046);
t2207 = sin(t2045) + sin(t2046);
t2038 = cos(t2055);
t2271 = 0.2e1 * t2038;
t1972 = t2018 * t2271 + t2011 * t2035 + (t2058 * t2107 + t2061 * t2108) * t2275 + (-t2107 * t2207 - t2108 * t2204) * pkin(3);
t2104 = cos(qJ(3,2));
t2265 = t2104 * pkin(3) + pkin(2);
t2020 = t2081 + t2265;
t2014 = 0.1e1 / t2020;
t2101 = sin(qJ(3,2));
t2111 = 0.2e1 * qJ(3,2);
t2086 = sin(t2111);
t2075 = pkin(3) * t2086;
t2005 = t2101 * t2287 + t2075 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1);
t1995 = 0.1e1 / t2005;
t2243 = t1972 * t1995;
t2176 = t2114 * t2243;
t2143 = t2069 * t2176;
t2002 = -t2035 * t2108 + t2038 * t2107;
t2233 = t2002 * t2014;
t2169 = t2101 * t2233;
t2211 = t2101 * t2104;
t2216 = t2070 * t2104;
t2015 = 0.1e1 / t2020 ^ 2;
t2232 = t2002 * t2015;
t2092 = t2104 ^ 2;
t2281 = MDP(6) * (0.2e1 * t2092 - 0.1e1);
t2289 = 0.2e1 * ((MDP(5) * t2211 + t2281) * t2232 - (MDP(10) * (t2104 * t2143 / 0.2e1 + t2070 * t2169) + MDP(11) * (t2216 * t2233 - t2101 * t2143 / 0.2e1)) * t2014) * t1972;
t2073 = qJ(1,1) + legFrame(1,3);
t2056 = pkin(7) + t2073;
t2036 = sin(t2056);
t2059 = sin(t2073);
t2062 = cos(t2073);
t2049 = qJ(3,1) + t2056;
t2050 = -qJ(3,1) + t2056;
t2203 = cos(t2049) + cos(t2050);
t2206 = sin(t2049) + sin(t2050);
t2039 = cos(t2056);
t2270 = 0.2e1 * t2039;
t1973 = t2018 * t2270 + t2011 * t2036 + (t2059 * t2107 + t2062 * t2108) * t2275 + (-t2107 * t2206 - t2108 * t2203) * pkin(3);
t2105 = cos(qJ(3,1));
t2264 = t2105 * pkin(3) + pkin(2);
t2021 = t2081 + t2264;
t2016 = 0.1e1 / t2021;
t2102 = sin(qJ(3,1));
t2112 = 0.2e1 * qJ(3,1);
t2087 = sin(t2112);
t2076 = pkin(3) * t2087;
t2006 = t2102 * t2287 + t2076 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1);
t1998 = 0.1e1 / t2006;
t2242 = t1973 * t1998;
t2174 = t2114 * t2242;
t2140 = t2069 * t2174;
t2003 = -t2036 * t2108 + t2039 * t2107;
t2231 = t2003 * t2016;
t2168 = t2102 * t2231;
t2209 = t2102 * t2105;
t2215 = t2070 * t2105;
t2017 = 0.1e1 / t2021 ^ 2;
t2230 = t2003 * t2017;
t2093 = t2105 ^ 2;
t2280 = MDP(6) * (0.2e1 * t2093 - 0.1e1);
t2288 = 0.2e1 * ((MDP(5) * t2209 + t2280) * t2230 - (MDP(10) * (t2105 * t2140 / 0.2e1 + t2070 * t2168) + MDP(11) * (t2215 * t2231 - t2102 * t2140 / 0.2e1)) * t2016) * t1973;
t2286 = MDP(4) / 0.2e1;
t2285 = MDP(10) / 0.2e1;
t2284 = MDP(11) / 0.2e1;
t2115 = 0.1e1 / pkin(3) ^ 2;
t2283 = t2115 / 0.2e1;
t2117 = pkin(1) ^ 2;
t2279 = MDP(4) * t2117 + MDP(1);
t2274 = -2 * pkin(2);
t1993 = 0.1e1 / t2004 ^ 2;
t1996 = 0.1e1 / t2005 ^ 2;
t1999 = 0.1e1 / t2006 ^ 2;
t2269 = -2 * t2109;
t2268 = 2 * t2109;
t2053 = pkin(6) + t2069;
t2267 = -t2053 / 0.2e1;
t2263 = pkin(3) * t2287;
t2261 = 2 * pkin(1);
t2113 = pkin(3) ^ 2;
t2135 = (pkin(2) ^ 2) + t2080 * t2268 + (t2109 ^ 2) + t2117;
t2202 = 0.2e1 * t2081;
t2007 = pkin(2) * t2202 + t2113 / 0.2e1 + t2135;
t2214 = t2100 * t2070;
t2008 = 0.1e1 / (t2074 + 0.2e1 * t2214);
t2145 = t2019 * t2178;
t2238 = t1992 * t2100;
t2179 = t1971 * t2238;
t2147 = t2053 * t2179;
t1947 = (-(0.2e1 * (-t2007 * t2235 + t2147) * t2103 + (-cos(t2110) + (-0.2e1 * t2091 - 0.1e1) * t2070 * t2012) * t2001 * pkin(3)) * t2235 + 0.2e1 * (t2085 * t2235 * t2267 + t2145) * t2244) * t2008;
t2259 = t1947 * t1992;
t2212 = t2101 * t2070;
t2009 = 0.1e1 / (t2075 + 0.2e1 * t2212);
t2142 = t2020 * t2176;
t2237 = t1995 * t2101;
t2177 = t1972 * t2237;
t2144 = t2053 * t2177;
t1948 = (-(0.2e1 * (-t2007 * t2233 + t2144) * t2104 + (-cos(t2111) + (-0.2e1 * t2092 - 0.1e1) * t2070 * t2014) * t2002 * pkin(3)) * t2233 + 0.2e1 * (t2086 * t2233 * t2267 + t2142) * t2243) * t2009;
t2258 = t1948 * t1995;
t2210 = t2102 * t2070;
t2010 = 0.1e1 / (t2076 + 0.2e1 * t2210);
t2139 = t2021 * t2174;
t2236 = t1998 * t2102;
t2175 = t1973 * t2236;
t2141 = t2053 * t2175;
t1949 = (-(0.2e1 * (-t2007 * t2231 + t2141) * t2105 + (-cos(t2112) + (-0.2e1 * t2093 - 0.1e1) * t2070 * t2016) * t2003 * pkin(3)) * t2231 + 0.2e1 * (t2087 * t2231 * t2267 + t2139) * t2242) * t2010;
t2257 = t1949 * t1998;
t1950 = -(-t2147 + (t2091 * t2113 + t2103 * t2263 + t2202 * t2266 + t2135) * t2235) / (t2074 / 0.2e1 + t2214) * t2114 * t2235 - 0.2e1 * (-t2053 * t2170 + t2103 * t2145) * t2008 * t2178;
t2256 = t1950 * t2103;
t1951 = -(-t2144 + (t2092 * t2113 + t2104 * t2263 + t2202 * t2265 + t2135) * t2233) / (t2075 / 0.2e1 + t2212) * t2114 * t2233 - 0.2e1 * (-t2053 * t2169 + t2104 * t2142) * t2009 * t2176;
t2255 = t1951 * t2104;
t1952 = -(-t2141 + (t2093 * t2113 + t2105 * t2263 + t2202 * t2264 + t2135) * t2231) / (t2076 / 0.2e1 + t2210) * t2114 * t2231 - 0.2e1 * (-t2053 * t2168 + t2105 * t2139) * t2010 * t2174;
t2254 = t1952 * t2105;
t2158 = t2109 + t2080;
t1953 = (-t2158 * t2235 + 0.2e1 * t2179) * t2234;
t2253 = t1953 * t2012;
t1954 = (-t2158 * t2233 + 0.2e1 * t2177) * t2232;
t2252 = t1954 * t2014;
t1955 = (-t2158 * t2231 + 0.2e1 * t2175) * t2230;
t2251 = t1955 * t2016;
t1968 = t1971 ^ 2;
t2250 = t1968 * t2100;
t2249 = t1968 * t2103;
t1969 = t1972 ^ 2;
t2248 = t1969 * t2101;
t2247 = t1969 * t2104;
t1970 = t1973 ^ 2;
t2246 = t1970 * t2102;
t2245 = t1970 * t2105;
t2241 = t2001 ^ 2 * t2013;
t2240 = t2002 ^ 2 * t2015;
t2239 = t2003 ^ 2 * t2017;
t2229 = t2012 * t2034;
t2228 = t2012 * t2037;
t2227 = t2014 * t2035;
t2226 = t2014 * t2038;
t2225 = t2016 * t2036;
t2224 = t2016 * t2039;
t2223 = t2069 * t2100;
t2222 = t2069 * t2101;
t2221 = t2069 * t2102;
t2220 = t2069 * t2103;
t2219 = t2069 * t2104;
t2218 = t2069 * t2105;
t2201 = 0.2e1 * MDP(6);
t2197 = t1950 * t2238;
t2196 = t1992 * t2256;
t2195 = t1951 * t2237;
t2194 = t1995 * t2255;
t2193 = t1952 * t2236;
t2192 = t1998 * t2254;
t2191 = t1993 * t2250;
t2190 = t1993 * t2249;
t1994 = t1992 * t1993;
t2189 = t1994 * t2250;
t2188 = t1994 * t2249;
t2187 = t1996 * t2248;
t2186 = t1996 * t2247;
t1997 = t1995 * t1996;
t2185 = t1997 * t2248;
t2184 = t1997 * t2247;
t2183 = t1999 * t2246;
t2182 = t1999 * t2245;
t2000 = t1998 * t1999;
t2181 = t2000 * t2246;
t2180 = t2000 * t2245;
t2173 = t2103 * t2241;
t2172 = t2104 * t2240;
t2171 = t2105 * t2239;
t2167 = t2100 ^ 2 * t2253;
t2166 = t2101 ^ 2 * t2252;
t2165 = t2102 ^ 2 * t2251;
t2164 = t1950 * t2229;
t2163 = t1951 * t2227;
t2162 = t1952 * t2225;
t2161 = t1950 * t2228;
t2160 = t1951 * t2226;
t2159 = t1952 * t2224;
t2153 = t2012 * t2191;
t2152 = t2012 * t2190;
t2151 = t2014 * t2187;
t2150 = t2014 * t2186;
t2149 = t2016 * t2183;
t2148 = t2016 * t2182;
t2138 = t2213 * t2253;
t2137 = t2211 * t2252;
t2136 = t2209 * t2251;
t2126 = -t2100 * MDP(5) * t2173 - t2241 * t2282 + t1950 * MDP(9) + (t1947 * t2103 - t1953 * t2223 + t2214 * t2241) * MDP(10) + (-t1947 * t2100 - t1953 * t2220 + t2070 * t2173) * MDP(11) + (MDP(7) * t2100 + MDP(8) * t2103) * t1953;
t2125 = -t2101 * MDP(5) * t2172 - t2240 * t2281 + t1951 * MDP(9) + (t1948 * t2104 - t1954 * t2222 + t2212 * t2240) * MDP(10) + (-t1948 * t2101 - t1954 * t2219 + t2070 * t2172) * MDP(11) + (MDP(7) * t2101 + MDP(8) * t2104) * t1954;
t2124 = -t2102 * MDP(5) * t2171 - t2239 * t2280 + t1952 * MDP(9) + (t1949 * t2105 - t1955 * t2221 + t2210 * t2239) * MDP(10) + (-t1949 * t2102 - t1955 * t2218 + t2070 * t2171) * MDP(11) + (MDP(7) * t2102 + MDP(8) * t2105) * t1955;
t2068 = -qJ(3,1) + t2073;
t2067 = qJ(3,1) + t2073;
t2066 = -qJ(3,2) + t2072;
t2065 = qJ(3,2) + t2072;
t2064 = -qJ(3,3) + t2071;
t2063 = qJ(3,3) + t2071;
t2051 = -0.2e1 * qJ(3,1) + t2056;
t2048 = t2112 + t2056;
t2047 = -0.2e1 * qJ(3,2) + t2055;
t2044 = t2111 + t2055;
t2043 = -0.2e1 * qJ(3,3) + t2054;
t2040 = t2110 + t2054;
t1979 = t2206 * t2287 + (sin(t2068) + sin(t2067)) * t2261 + t2203 * t2269 + (sin(t2051) + sin(t2048) + 0.2e1 * t2036) * pkin(3);
t1978 = t2207 * t2287 + (sin(t2066) + sin(t2065)) * t2261 + t2204 * t2269 + (sin(t2047) + sin(t2044) + 0.2e1 * t2035) * pkin(3);
t1977 = t2208 * t2287 + (sin(t2064) + sin(t2063)) * t2261 + t2205 * t2269 + (sin(t2043) + sin(t2040) + 0.2e1 * t2034) * pkin(3);
t1976 = t2203 * t2287 + (cos(t2068) + cos(t2067)) * t2261 + t2206 * t2268 + (cos(t2051) + cos(t2048) + t2270) * pkin(3);
t1975 = t2204 * t2287 + (cos(t2066) + cos(t2065)) * t2261 + t2207 * t2268 + (cos(t2047) + cos(t2044) + t2271) * pkin(3);
t1974 = t2205 * t2287 + (cos(t2064) + cos(t2063)) * t2261 + t2208 * t2268 + (cos(t2043) + cos(t2040) + t2272) * pkin(3);
t1946 = -t1952 * t2221 + 0.2e1 * t1955 * t2215;
t1945 = -t1950 * t2223 + 0.2e1 * t1953 * t2217;
t1944 = -t1951 * t2222 + 0.2e1 * t1954 * t2216;
t1943 = -t1952 * t2218 - 0.2e1 * t1955 * t2210;
t1942 = -t1951 * t2219 - 0.2e1 * t1954 * t2212;
t1941 = -t1950 * t2220 - 0.2e1 * t1953 * t2214;
t1 = [(-t2034 * t2167 - t2035 * t2166 - t2036 * t2165) * MDP(5) + (-t2034 * t2138 - t2035 * t2137 - t2036 * t2136) * t2201 + (-t2100 * t2164 - t2101 * t2163 - t2102 * t2162) * MDP(7) + (-t2103 * t2164 - t2104 * t2163 - t2105 * t2162) * MDP(8) + (-t1944 * t2227 - t1945 * t2229 - t1946 * t2225) * MDP(10) + (-t1941 * t2229 - t1942 * t2227 - t1943 * t2225) * MDP(11) + ((-t2034 * t2152 - t2035 * t2150 - t2036 * t2148) * MDP(7) + (t2034 * t2153 + t2035 * t2151 + t2036 * t2149) * MDP(8)) * t2115 + (t1974 * t2259 + t1975 * t2258 + t1976 * t2257) * t2286 + (t1974 * t2196 + t1975 * t2194 + t1976 * t2192) * t2285 + (-t1974 * t2197 - t1975 * t2195 - t1976 * t2193) * t2284 + ((-t1974 * t2189 - t1975 * t2185 - t1976 * t2181) * MDP(10) + (-t1974 * t2188 - t1975 * t2184 - t1976 * t2180) * MDP(11)) * t2283 + ((-t2036 * t2288 + t2124 * (-pkin(3) * t2203 + t2036 * t2269 + t2039 * t2274 + t2062 * t2275)) * t1998 + (-t2035 * t2289 + t2125 * (-pkin(3) * t2204 + t2035 * t2269 + t2038 * t2274 + t2061 * t2275)) * t1995 + (-t2034 * t2290 + t2126 * (-pkin(3) * t2205 + t2034 * t2269 + t2037 * t2274 + t2060 * t2275)) * t1992) * t2114 + t2279 * (-t1953 * t2229 - t1954 * t2227 - t1955 * t2225); (t2037 * t2167 + t2038 * t2166 + t2039 * t2165) * MDP(5) + (t2037 * t2138 + t2038 * t2137 + t2039 * t2136) * t2201 + (t2100 * t2161 + t2101 * t2160 + t2102 * t2159) * MDP(7) + (t2103 * t2161 + t2104 * t2160 + t2105 * t2159) * MDP(8) + (t1944 * t2226 + t1945 * t2228 + t1946 * t2224) * MDP(10) + (t1941 * t2228 + t1942 * t2226 + t1943 * t2224) * MDP(11) + ((t2037 * t2152 + t2038 * t2150 + t2039 * t2148) * MDP(7) + (-t2037 * t2153 - t2038 * t2151 - t2039 * t2149) * MDP(8)) * t2115 + (t1977 * t2259 + t1978 * t2258 + t1979 * t2257) * t2286 + (t1977 * t2196 + t1978 * t2194 + t1979 * t2192) * t2285 + (-t1977 * t2197 - t1978 * t2195 - t1979 * t2193) * t2284 + ((-t1977 * t2189 - t1978 * t2185 - t1979 * t2181) * MDP(10) + (-t1977 * t2188 - t1978 * t2184 - t1979 * t2180) * MDP(11)) * t2283 + ((t2039 * t2288 + t2124 * (-pkin(3) * t2206 + t2036 * t2274 + t2039 * t2268 + t2059 * t2275)) * t1998 + (t2038 * t2289 + t2125 * (-pkin(3) * t2207 + t2035 * t2274 + t2038 * t2268 + t2058 * t2275)) * t1995 + (t2037 * t2290 + t2126 * (-pkin(3) * t2208 + t2034 * t2274 + t2037 * t2268 + t2057 * t2275)) * t1992) * t2114 + t2279 * (t1953 * t2228 + t1954 * t2226 + t1955 * t2224); (t1947 + t1948 + t1949) * MDP(4) + (t2254 + t2255 + t2256) * MDP(10) + (-t1950 * t2100 - t1951 * t2101 - t1952 * t2102) * MDP(11) + ((-t2183 - t2187 - t2191) * MDP(10) + (-t2182 - t2186 - t2190) * MDP(11)) * t2115;];
taucX  = t1;
