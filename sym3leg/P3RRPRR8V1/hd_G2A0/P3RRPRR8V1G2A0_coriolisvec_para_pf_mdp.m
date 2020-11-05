% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:50
% EndTime: 2020-08-06 19:59:55
% DurationCPUTime: 5.66s
% Computational Cost: add. (19026->334), mult. (42941->669), div. (2031->24), fcn. (28644->35), ass. (0->296)
t2099 = pkin(4) + qJ(3,3);
t2082 = 0.1e1 / t2099;
t2135 = (t2099 ^ 2);
t2083 = 0.1e1 / t2135;
t2084 = t2082 * t2083;
t2100 = pkin(4) + qJ(3,2);
t2085 = 0.1e1 / t2100;
t2137 = (t2100 ^ 2);
t2086 = 0.1e1 / t2137;
t2087 = t2085 * t2086;
t2101 = pkin(4) + qJ(3,1);
t2088 = 0.1e1 / t2101;
t2139 = (t2101 ^ 2);
t2089 = 0.1e1 / t2139;
t2090 = t2088 * t2089;
t2308 = 2 * pkin(1);
t2097 = sin(pkin(5));
t2119 = xDP(2);
t2297 = pkin(2) * t2119;
t2051 = t2097 * t2297;
t2120 = xDP(1);
t2296 = pkin(2) * t2120;
t2052 = t2097 * t2296;
t2106 = sin(qJ(1,3));
t2056 = t2106 * t2099;
t2102 = legFrame(3,2);
t2068 = sin(t2102);
t2071 = cos(t2102);
t2105 = sin(qJ(2,3));
t2111 = cos(qJ(2,3));
t2112 = cos(qJ(1,3));
t2118 = xDP(3);
t2098 = cos(pkin(5));
t2298 = pkin(2) * t2098;
t2055 = pkin(1) + t2298;
t2171 = t2055 * t2120;
t2299 = pkin(2) * t2097;
t2220 = t2118 * t2299;
t2183 = t2112 * t2220;
t2188 = t2106 * t2052;
t2189 = t2106 * t2051;
t2239 = t2099 * t2112;
t2246 = t2055 * t2119;
t2250 = t2055 * t2111;
t1987 = ((t2106 * t2171 + t2051) * t2111 - t2120 * t2239) * t2071 + ((-t2106 * t2246 + t2052) * t2111 + t2119 * t2239) * t2068 + t2118 * (t2112 * t2250 + t2056) + ((-t2188 + t2246) * t2071 + (t2189 + t2171) * t2068 - t2183) * t2105;
t2320 = 0.2e1 * t1987;
t2108 = sin(qJ(1,2));
t2057 = t2108 * t2100;
t2103 = legFrame(2,2);
t2069 = sin(t2103);
t2072 = cos(t2103);
t2107 = sin(qJ(2,2));
t2113 = cos(qJ(2,2));
t2114 = cos(qJ(1,2));
t2182 = t2114 * t2220;
t2186 = t2108 * t2052;
t2187 = t2108 * t2051;
t2238 = t2100 * t2114;
t2249 = t2055 * t2113;
t1988 = ((t2108 * t2171 + t2051) * t2113 - t2120 * t2238) * t2072 + ((-t2108 * t2246 + t2052) * t2113 + t2119 * t2238) * t2069 + t2118 * (t2114 * t2249 + t2057) + ((-t2186 + t2246) * t2072 + (t2187 + t2171) * t2069 - t2182) * t2107;
t2319 = 0.2e1 * t1988;
t2110 = sin(qJ(1,1));
t2058 = t2110 * t2101;
t2104 = legFrame(1,2);
t2070 = sin(t2104);
t2073 = cos(t2104);
t2109 = sin(qJ(2,1));
t2115 = cos(qJ(2,1));
t2116 = cos(qJ(1,1));
t2181 = t2116 * t2220;
t2184 = t2110 * t2052;
t2185 = t2110 * t2051;
t2237 = t2101 * t2116;
t2248 = t2055 * t2115;
t1989 = ((t2110 * t2171 + t2051) * t2115 - t2120 * t2237) * t2073 + ((-t2110 * t2246 + t2052) * t2115 + t2119 * t2237) * t2070 + t2118 * (t2116 * t2248 + t2058) + ((-t2184 + t2246) * t2073 + (t2185 + t2171) * t2070 - t2181) * t2109;
t2318 = 0.2e1 * t1989;
t2029 = t2068 * t2120 + t2071 * t2119;
t2079 = qJ(2,3) + pkin(5);
t2044 = pkin(1) * t2105 + pkin(2) * sin(t2079);
t2124 = pkin(2) ^ 2;
t2125 = pkin(1) ^ 2;
t2225 = -t2124 - t2125;
t2050 = t2298 * t2308 - t2225;
t2121 = 0.2e1 * qJ(2,3);
t2224 = pkin(2) * t2308;
t2074 = t2111 * pkin(1);
t2314 = t2074 + pkin(2) * cos(t2079);
t2032 = 0.1e1 / t2314;
t2268 = t2029 * t2032;
t2226 = pkin(1) * t2120 + t2098 * t2296;
t2227 = pkin(1) * t2119 + t2098 * t2297;
t2247 = t2055 * t2118;
t1993 = (t2112 * t2247 + t2051 * t2071 + t2052 * t2068 + (-t2068 * t2246 + t2071 * t2171) * t2106) * t2111 + t2105 * ((-t2188 + t2227) * t2071 + (t2189 + t2226) * t2068 - t2183);
t2242 = t2097 * t2105;
t2023 = -pkin(2) * t2242 + t2250;
t2014 = 0.1e1 / t2023;
t2280 = t1993 * t2014;
t2286 = t1987 * t2084;
t2146 = -(t2044 * t2280 - t2050 * t2268) * t2082 * t2268 + (-t2044 * t2099 * t2268 * t2083 - t2314 * t2286 - ((-t2124 * cos(0.2e1 * t2079) - cos(t2121) * t2125 - (2 * t2135) + (-cos(t2121 + pkin(5)) - t2098) * t2224 + t2225) * t2280 + t2314 * t2320) * t2084 / 0.2e1) * t2280;
t2206 = t2105 * t2280;
t2180 = t2032 * t2206;
t2166 = t2082 * t2180;
t2015 = 0.1e1 / t2023 ^ 2;
t2283 = t1993 ^ 2 * t2015;
t2213 = t2083 * t2283;
t2026 = t2029 ^ 2;
t2126 = t2314 ^ 2;
t2033 = 0.1e1 / t2126;
t2205 = t2111 * t2280;
t2279 = t1993 * t2083;
t2207 = t2014 * t2279;
t1966 = -t2280 * t2286 + (-(-t2055 * t2205 + t2206 * t2299 + t1987) * t2207 + t2026 * t2033 * t2050 / (t2074 + (t2098 * t2111 - t2242) * pkin(2))) * t2082;
t2293 = t1966 * t2111;
t2317 = t2029 * t2166 * t2308 - pkin(1) * t2293 - qJ(3,3) * t2213 + t2146;
t2030 = t2069 * t2120 + t2072 * t2119;
t2080 = qJ(2,2) + pkin(5);
t2045 = pkin(1) * t2107 + pkin(2) * sin(t2080);
t2122 = 0.2e1 * qJ(2,2);
t2075 = t2113 * pkin(1);
t2313 = t2075 + pkin(2) * cos(t2080);
t2036 = 0.1e1 / t2313;
t2267 = t2030 * t2036;
t1994 = (t2114 * t2247 + t2051 * t2072 + t2052 * t2069 + (-t2069 * t2246 + t2072 * t2171) * t2108) * t2113 + t2107 * ((-t2186 + t2227) * t2072 + (t2187 + t2226) * t2069 - t2182);
t2241 = t2097 * t2107;
t2024 = -pkin(2) * t2241 + t2249;
t2016 = 0.1e1 / t2024;
t2278 = t1994 * t2016;
t2285 = t1988 * t2087;
t2145 = -(t2045 * t2278 - t2050 * t2267) * t2085 * t2267 + (-t2045 * t2100 * t2267 * t2086 - t2313 * t2285 - ((-t2124 * cos(0.2e1 * t2080) - cos(t2122) * t2125 - (2 * t2137) + (-cos(pkin(5) + t2122) - t2098) * t2224 + t2225) * t2278 + t2313 * t2319) * t2087 / 0.2e1) * t2278;
t2203 = t2107 * t2278;
t2178 = t2036 * t2203;
t2165 = t2085 * t2178;
t2017 = 0.1e1 / t2024 ^ 2;
t2282 = t1994 ^ 2 * t2017;
t2211 = t2086 * t2282;
t2027 = t2030 ^ 2;
t2129 = t2313 ^ 2;
t2037 = 0.1e1 / t2129;
t2202 = t2113 * t2278;
t2277 = t1994 * t2086;
t2204 = t2016 * t2277;
t1967 = -t2278 * t2285 + (-(-t2055 * t2202 + t2203 * t2299 + t1988) * t2204 + t2027 * t2037 * t2050 / (t2075 + (t2098 * t2113 - t2241) * pkin(2))) * t2085;
t2290 = t1967 * t2113;
t2316 = t2030 * t2165 * t2308 - pkin(1) * t2290 - qJ(3,2) * t2211 + t2145;
t2031 = t2070 * t2120 + t2073 * t2119;
t2081 = qJ(2,1) + pkin(5);
t2046 = pkin(1) * t2109 + pkin(2) * sin(t2081);
t2123 = 0.2e1 * qJ(2,1);
t2076 = t2115 * pkin(1);
t2312 = t2076 + pkin(2) * cos(t2081);
t2040 = 0.1e1 / t2312;
t2266 = t2031 * t2040;
t1995 = (t2116 * t2247 + t2051 * t2073 + t2052 * t2070 + (-t2070 * t2246 + t2073 * t2171) * t2110) * t2115 + t2109 * ((-t2184 + t2227) * t2073 + (t2185 + t2226) * t2070 - t2181);
t2240 = t2097 * t2109;
t2025 = -pkin(2) * t2240 + t2248;
t2018 = 0.1e1 / t2025;
t2276 = t1995 * t2018;
t2284 = t1989 * t2090;
t2144 = -(t2046 * t2276 - t2050 * t2266) * t2088 * t2266 + (-t2046 * t2101 * t2266 * t2089 - t2312 * t2284 - ((-t2124 * cos(0.2e1 * t2081) - cos(t2123) * t2125 - (2 * t2139) + (-cos(pkin(5) + t2123) - t2098) * t2224 + t2225) * t2276 + t2312 * t2318) * t2090 / 0.2e1) * t2276;
t2200 = t2109 * t2276;
t2176 = t2040 * t2200;
t2164 = t2088 * t2176;
t2019 = 0.1e1 / t2025 ^ 2;
t2281 = t1995 ^ 2 * t2019;
t2209 = t2089 * t2281;
t2028 = t2031 ^ 2;
t2132 = t2312 ^ 2;
t2041 = 0.1e1 / t2132;
t2199 = t2115 * t2276;
t2275 = t1995 * t2089;
t2201 = t2018 * t2275;
t1968 = -t2276 * t2284 + (-(-t2055 * t2199 + t2200 * t2299 + t1989) * t2201 + t2028 * t2041 * t2050 / (t2076 + (t2098 * t2115 - t2240) * pkin(2))) * t2088;
t2287 = t1968 * t2115;
t2315 = t2031 * t2164 * t2308 - pkin(1) * t2287 - qJ(3,1) * t2209 + t2144;
t2311 = t2317 * t2082;
t2310 = t2316 * t2085;
t2309 = t2315 * t2088;
t2307 = 2 * MDP(5);
t2306 = -2 * t2125;
t2094 = t2111 ^ 2;
t2065 = 0.2e1 * t2094 - 0.1e1;
t2095 = t2113 ^ 2;
t2066 = 0.2e1 * t2095 - 0.1e1;
t2096 = t2115 ^ 2;
t2067 = 0.2e1 * t2096 - 0.1e1;
t2305 = pkin(1) * t2026;
t2304 = pkin(1) * t2027;
t2303 = pkin(1) * t2028;
t2295 = t1966 * t2082;
t2294 = t1966 * t2105;
t2292 = t1967 * t2085;
t2291 = t1967 * t2107;
t2289 = t1968 * t2088;
t2288 = t1968 * t2109;
t2265 = t2032 * t2068;
t2264 = t2032 * t2071;
t2263 = t2033 * t2111;
t2262 = t2032 * t2033 * t2044;
t2261 = t2036 * t2069;
t2260 = t2036 * t2072;
t2259 = t2037 * t2113;
t2258 = t2036 * t2037 * t2045;
t2257 = t2040 * t2070;
t2256 = t2040 * t2073;
t2255 = t2041 * t2115;
t2254 = t2040 * t2041 * t2046;
t2253 = t2055 * t2106;
t2252 = t2055 * t2108;
t2251 = t2055 * t2110;
t2245 = t2082 * t2112;
t2244 = t2085 * t2114;
t2243 = t2088 * t2116;
t2236 = t2105 * t2111;
t2235 = t2107 * t2113;
t2234 = t2109 * t2115;
t1978 = t2207 * t2320 - t2263 * t2305;
t2157 = t2029 * t2111 * t2180;
t2192 = t2105 * t2262;
t2174 = t2026 * t2192;
t2233 = (qJ(3,3) ^ 2 + t2094 * t2125) * t1966 + (-qJ(3,3) * t2174 - t2111 * t2146) * pkin(1) + t2082 * t2157 * t2306 + t1978 * qJ(3,3);
t1979 = t2204 * t2319 - t2259 * t2304;
t2156 = t2030 * t2113 * t2178;
t2191 = t2107 * t2258;
t2173 = t2027 * t2191;
t2232 = (qJ(3,2) ^ 2 + t2095 * t2125) * t1967 + (-qJ(3,2) * t2173 - t2113 * t2145) * pkin(1) + t2085 * t2156 * t2306 + t1979 * qJ(3,2);
t1980 = t2201 * t2318 - t2255 * t2303;
t2155 = t2031 * t2115 * t2176;
t2190 = t2109 * t2254;
t2172 = t2028 * t2190;
t2231 = (qJ(3,1) ^ 2 + t2096 * t2125) * t1968 + (-qJ(3,1) * t2172 - t2115 * t2144) * pkin(1) + t2088 * t2155 * t2306 + t1980 * qJ(3,1);
t2230 = -pkin(1) * t2174 + 0.2e1 * qJ(3,3) * t1966 + t1978;
t2229 = -pkin(1) * t2173 + 0.2e1 * qJ(3,2) * t1967 + t1979;
t2228 = -pkin(1) * t2172 + 0.2e1 * qJ(3,1) * t1968 + t1980;
t2223 = t2106 * t2299;
t2222 = t2108 * t2299;
t2221 = t2110 * t2299;
t2219 = t2032 * t2294;
t2218 = t2032 * t2293;
t2217 = t2036 * t2291;
t2216 = t2036 * t2290;
t2215 = t2040 * t2288;
t2214 = t2040 * t2287;
t2212 = t2084 * t2283;
t2210 = t2087 * t2282;
t2208 = t2090 * t2281;
t2198 = t2026 / t2126 ^ 2 * t2044;
t2197 = t2026 * t2245;
t2196 = t2027 / t2129 ^ 2 * t2045;
t2195 = t2027 * t2244;
t2194 = t2028 / t2132 ^ 2 * t2046;
t2193 = t2028 * t2243;
t2179 = t2268 * t2279;
t2177 = t2267 * t2277;
t2175 = t2266 * t2275;
t2163 = (pkin(1) * t2205 - 0.2e1 * t1987) * t2082 * t2166;
t2162 = (pkin(1) * t2202 - 0.2e1 * t1988) * t2085 * t2165;
t2161 = (pkin(1) * t2199 - 0.2e1 * t1989) * t2088 * t2164;
t2160 = t2032 * t2213 * t2236;
t2159 = t2036 * t2211 * t2235;
t2158 = t2040 * t2209 * t2234;
t2154 = -t2033 * t2105 + t2111 * t2262;
t2153 = t2192 + t2263;
t2152 = -t2037 * t2107 + t2113 * t2258;
t2151 = t2191 + t2259;
t2150 = -t2041 * t2109 + t2115 * t2254;
t2149 = t2190 + t2255;
t2148 = t2068 * t2219 + t2069 * t2217 + t2070 * t2215;
t2147 = t2071 * t2219 + t2072 * t2217 + t2073 * t2215;
t2091 = t2105 ^ 2;
t2143 = 0.2e1 * (MDP(4) * t2236 + MDP(5) * t2065) * t2015 * t2179 + (t2230 * MDP(11) + t2233 * MDP(12) + (t2091 * MDP(4) + t2236 * t2307 + MDP(1)) * t1966 + (t2153 * MDP(6) + t2154 * MDP(7)) * t2026) * t2082 * t2014;
t2092 = t2107 ^ 2;
t2142 = 0.2e1 * (MDP(4) * t2235 + MDP(5) * t2066) * t2017 * t2177 + (t2229 * MDP(11) + t2232 * MDP(12) + (t2092 * MDP(4) + t2235 * t2307 + MDP(1)) * t1967 + (t2151 * MDP(6) + t2152 * MDP(7)) * t2027) * t2085 * t2016;
t2093 = t2109 ^ 2;
t2141 = 0.2e1 * (MDP(4) * t2234 + MDP(5) * t2067) * t2019 * t2175 + (t2228 * MDP(11) + t2231 * MDP(12) + (t2093 * MDP(4) + t2234 * t2307 + MDP(1)) * t1968 + (t2149 * MDP(6) + t2150 * MDP(7)) * t2028) * t2088 * t2018;
t2022 = t2055 * t2109 + t2115 * t2299;
t2021 = t2055 * t2107 + t2113 * t2299;
t2020 = t2055 * t2105 + t2111 * t2299;
t2013 = t2116 * t2312 + t2058;
t2012 = t2114 * t2313 + t2057;
t2011 = t2112 * t2314 + t2056;
t2010 = t2025 * t2110 - t2237;
t2009 = t2024 * t2108 - t2238;
t2008 = t2023 * t2106 - t2239;
t2001 = -t2010 * t2070 + t2022 * t2073;
t2000 = t2010 * t2073 + t2022 * t2070;
t1999 = -t2009 * t2069 + t2021 * t2072;
t1998 = t2009 * t2072 + t2021 * t2069;
t1997 = -t2008 * t2068 + t2020 * t2071;
t1996 = t2008 * t2071 + t2020 * t2068;
t1983 = t2067 * t2209;
t1982 = t2066 * t2211;
t1981 = t2065 * t2213;
t1962 = -qJ(3,1) * t2288 + t2254 * t2303;
t1961 = -qJ(3,2) * t2291 + t2258 * t2304;
t1960 = -qJ(3,3) * t2294 + t2262 * t2305;
t1 = [(-t2068 * t2160 - t2069 * t2159 - t2070 * t2158) * MDP(4) + (-t1981 * t2265 - t1982 * t2261 - t1983 * t2257) * MDP(5) + t2148 * MDP(6) + (t2068 * t2218 + t2069 * t2216 + t2070 * t2214) * MDP(7) + (t2068 * t2198 + t2069 * t2196 + t2070 * t2194) * MDP(8) + (-t1996 * t2212 - t1998 * t2210 - t2000 * t2208) * MDP(11) + (t1996 * t2311 + t1998 * t2310 + t2000 * t2309) * MDP(12) + t2141 * ((t2070 * t2299 + t2073 * t2251) * t2115 + (t2055 * t2070 - t2073 * t2221) * t2109) + t2142 * ((t2069 * t2299 + t2072 * t2252) * t2113 + (t2055 * t2069 - t2072 * t2222) * t2107) + t2143 * ((t2068 * t2299 + t2071 * t2253) * t2111 + (t2055 * t2068 - t2071 * t2223) * t2105) + (-t2148 * MDP(11) + (t1960 * t2265 + t1961 * t2261 + t1962 * t2257 + t2068 * t2163 + t2069 * t2162 + t2070 * t2161) * MDP(12)) * pkin(1); (-t2071 * t2160 - t2072 * t2159 - t2073 * t2158) * MDP(4) + (-t1981 * t2264 - t1982 * t2260 - t1983 * t2256) * MDP(5) + t2147 * MDP(6) + (t2071 * t2218 + t2072 * t2216 + t2073 * t2214) * MDP(7) + (t2071 * t2198 + t2072 * t2196 + t2073 * t2194) * MDP(8) + (-t1997 * t2212 - t1999 * t2210 - t2001 * t2208) * MDP(11) + (t1997 * t2311 + t1999 * t2310 + t2001 * t2309) * MDP(12) + t2141 * ((-t2070 * t2251 + t2073 * t2299) * t2115 + t2109 * (t2055 * t2073 + t2070 * t2221)) + t2142 * ((-t2069 * t2252 + t2072 * t2299) * t2113 + t2107 * (t2055 * t2072 + t2069 * t2222)) + t2143 * ((-t2068 * t2253 + t2071 * t2299) * t2111 + t2105 * (t2055 * t2071 + t2068 * t2223)) + (-t2147 * MDP(11) + (t1960 * t2264 + t1961 * t2260 + t1962 * t2256 + t2071 * t2163 + t2072 * t2162 + t2073 * t2161) * MDP(12)) * pkin(1); (t1966 * t2245 + t1967 * t2244 + t1968 * t2243) * MDP(1) + ((0.2e1 * t2089 * t2155 + t2093 * t2289) * t2116 + (0.2e1 * t2086 * t2156 + t2092 * t2292) * t2114 + (0.2e1 * t2083 * t2157 + t2091 * t2295) * t2112) * MDP(4) + ((t2018 * t2067 * t2175 + t2234 * t2289) * t2116 + (t2016 * t2066 * t2177 + t2235 * t2292) * t2114 + (t2014 * t2065 * t2179 + t2236 * t2295) * t2112) * t2307 + (t2149 * t2193 + t2151 * t2195 + t2153 * t2197) * MDP(6) + (t2150 * t2193 + t2152 * t2195 + t2154 * t2197) * MDP(7) + (-t2011 * t2212 - t2012 * t2210 - t2013 * t2208 + t2228 * t2243 + t2229 * t2244 + t2230 * t2245) * MDP(11) + ((t2315 * t2013 + t2231 * t2116) * t2088 + (t2316 * t2012 + t2232 * t2114) * t2085 + (t2317 * t2011 + t2233 * t2112) * t2082) * MDP(12);];
taucX  = t1;