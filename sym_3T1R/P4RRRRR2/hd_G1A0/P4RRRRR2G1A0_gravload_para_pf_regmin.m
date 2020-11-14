% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4RRRRR2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x17]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4RRRRR2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:25:45
% EndTime: 2020-08-07 17:25:47
% DurationCPUTime: 1.62s
% Computational Cost: add. (2033->261), mult. (1756->432), div. (484->18), fcn. (1804->90), ass. (0->222)
t2180 = legFrame(4,3);
t2154 = sin(t2180);
t2287 = cos(t2180);
t2084 = t2154 * g(1) - t2287 * g(2);
t2088 = t2287 * g(1) + t2154 * g(2);
t2160 = qJ(1,4) + qJ(2,4);
t2146 = sin(t2160);
t2147 = cos(t2160);
t2056 = t2084 * t2147 + t2088 * t2146;
t2184 = sin(qJ(3,4));
t2292 = t2056 * t2184 ^ 2;
t2181 = legFrame(3,3);
t2155 = sin(t2181);
t2286 = cos(t2181);
t2085 = t2155 * g(1) - t2286 * g(2);
t2089 = t2286 * g(1) + t2155 * g(2);
t2177 = qJ(1,3) + qJ(2,3);
t2148 = sin(t2177);
t2151 = cos(t2177);
t2058 = t2085 * t2151 + t2089 * t2148;
t2188 = sin(qJ(3,3));
t2291 = t2058 * t2188 ^ 2;
t2182 = legFrame(2,3);
t2156 = sin(t2182);
t2285 = cos(t2182);
t2086 = t2156 * g(1) - t2285 * g(2);
t2090 = t2285 * g(1) + t2156 * g(2);
t2178 = qJ(1,2) + qJ(2,2);
t2149 = sin(t2178);
t2152 = cos(t2178);
t2059 = t2086 * t2152 + t2090 * t2149;
t2190 = sin(qJ(3,2));
t2290 = t2059 * t2190 ^ 2;
t2183 = legFrame(1,3);
t2157 = sin(t2183);
t2284 = cos(t2183);
t2087 = t2157 * g(1) - t2284 * g(2);
t2091 = t2284 * g(1) + t2157 * g(2);
t2179 = qJ(1,1) + qJ(2,1);
t2150 = sin(t2179);
t2153 = cos(t2179);
t2060 = t2087 * t2153 + t2091 * t2150;
t2192 = sin(qJ(3,1));
t2289 = t2060 * t2192 ^ 2;
t2288 = 2 * pkin(1);
t2243 = qJ(1,4) + t2180;
t2283 = pkin(1) * sin(t2243);
t2277 = qJ(1,3) + t2181;
t2282 = pkin(1) * sin(t2277);
t2278 = qJ(1,2) + t2182;
t2281 = pkin(1) * sin(t2278);
t2279 = qJ(1,1) + t2183;
t2280 = pkin(1) * sin(t2279);
t2276 = t2056 * t2184;
t2186 = cos(qJ(3,4));
t2275 = t2056 * t2186;
t2274 = t2058 * t2188;
t2194 = cos(qJ(3,3));
t2273 = t2058 * t2194;
t2272 = t2059 * t2190;
t2196 = cos(qJ(3,2));
t2271 = t2059 * t2196;
t2270 = t2060 * t2192;
t2198 = cos(qJ(3,1));
t2269 = t2060 * t2198;
t2135 = qJ(2,4) + t2243;
t2120 = qJ(3,4) + t2135;
t2099 = sin(t2120);
t2121 = -qJ(3,4) + t2135;
t2100 = sin(t2121);
t2072 = -0.2e1 * t2283 + (-t2099 - t2100) * pkin(2);
t2094 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t2268 = t2072 * t2094;
t2073 = cos(t2243) * t2288 + (cos(t2120) + cos(t2121)) * pkin(2);
t2267 = t2073 * t2094;
t2139 = qJ(2,3) + t2277;
t2128 = qJ(3,3) + t2139;
t2105 = sin(t2128);
t2129 = -qJ(3,3) + t2139;
t2106 = sin(t2129);
t2074 = -0.2e1 * t2282 + (-t2105 - t2106) * pkin(2);
t2095 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t2266 = t2074 * t2095;
t2140 = qJ(2,2) + t2278;
t2130 = qJ(3,2) + t2140;
t2107 = sin(t2130);
t2131 = -qJ(3,2) + t2140;
t2108 = sin(t2131);
t2075 = -0.2e1 * t2281 + (-t2107 - t2108) * pkin(2);
t2096 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t2265 = t2075 * t2096;
t2141 = qJ(2,1) + t2279;
t2132 = qJ(3,1) + t2141;
t2109 = sin(t2132);
t2133 = -qJ(3,1) + t2141;
t2110 = sin(t2133);
t2076 = -0.2e1 * t2280 + (-t2109 - t2110) * pkin(2);
t2097 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t2264 = t2076 * t2097;
t2077 = cos(t2277) * t2288 + (cos(t2128) + cos(t2129)) * pkin(2);
t2263 = t2077 * t2095;
t2078 = cos(t2278) * t2288 + (cos(t2130) + cos(t2131)) * pkin(2);
t2262 = t2078 * t2096;
t2079 = cos(t2279) * t2288 + (cos(t2132) + cos(t2133)) * pkin(2);
t2261 = t2079 * t2097;
t2200 = xP(4);
t2104 = -t2200 + t2135;
t2201 = koppelP(4,2);
t2205 = koppelP(4,1);
t2080 = -t2201 * cos(t2104) + t2205 * sin(t2104);
t2162 = 0.1e1 / sin(qJ(2,4));
t2260 = t2080 * t2162;
t2112 = -t2200 + t2139;
t2202 = koppelP(3,2);
t2206 = koppelP(3,1);
t2081 = -t2202 * cos(t2112) + t2206 * sin(t2112);
t2166 = 0.1e1 / sin(qJ(2,3));
t2259 = t2081 * t2166;
t2113 = -t2200 + t2140;
t2203 = koppelP(2,2);
t2207 = koppelP(2,1);
t2082 = -t2203 * cos(t2113) + t2207 * sin(t2113);
t2168 = 0.1e1 / sin(qJ(2,2));
t2258 = t2082 * t2168;
t2114 = -t2200 + t2141;
t2204 = koppelP(1,2);
t2208 = koppelP(1,1);
t2083 = -t2204 * cos(t2114) + t2208 * sin(t2114);
t2170 = 0.1e1 / sin(qJ(2,1));
t2257 = t2083 * t2170;
t2118 = sin(t2135);
t2256 = t2118 * t2162;
t2119 = cos(t2135);
t2255 = t2119 * t2162;
t2122 = sin(t2139);
t2254 = t2122 * t2166;
t2123 = sin(t2140);
t2253 = t2123 * t2168;
t2124 = sin(t2141);
t2252 = t2124 * t2170;
t2125 = cos(t2139);
t2251 = t2125 * t2166;
t2126 = cos(t2140);
t2250 = t2126 * t2168;
t2127 = cos(t2141);
t2249 = t2127 * t2170;
t2248 = t2162 * t2184;
t2247 = t2166 * t2188;
t2246 = t2168 * t2190;
t2245 = t2170 * t2192;
t2209 = 0.1e1 / pkin(2);
t2210 = 1 / pkin(1);
t2244 = t2209 * t2210;
t2242 = t2094 * t2276;
t2241 = t2094 * t2275;
t2240 = t2162 * t2275;
t2239 = t2095 * t2274;
t2238 = t2095 * t2273;
t2237 = t2166 * t2273;
t2236 = t2096 * t2272;
t2235 = t2096 * t2271;
t2234 = t2168 * t2271;
t2233 = t2097 * t2270;
t2232 = t2097 * t2269;
t2231 = t2170 * t2269;
t2098 = cos(qJ(2,4)) * pkin(1) + t2186 * pkin(2);
t2230 = t2098 * t2162 / t2186 ^ 2;
t2101 = cos(qJ(2,3)) * pkin(1) + t2194 * pkin(2);
t2229 = t2101 * t2166 / t2194 ^ 2;
t2102 = cos(qJ(2,2)) * pkin(1) + t2196 * pkin(2);
t2228 = t2102 * t2168 / t2196 ^ 2;
t2103 = cos(qJ(2,1)) * pkin(1) + t2198 * pkin(2);
t2227 = t2103 * t2170 / t2198 ^ 2;
t2163 = 0.1e1 / t2186;
t2226 = t2163 * t2248;
t2171 = 0.1e1 / t2194;
t2225 = t2171 * t2247;
t2173 = 0.1e1 / t2196;
t2224 = t2173 * t2246;
t2175 = 0.1e1 / t2198;
t2223 = t2175 * t2245;
t2222 = t2056 * t2248;
t2221 = t2058 * t2247;
t2220 = t2059 * t2246;
t2219 = t2060 * t2245;
t2218 = t2163 * t2222;
t2217 = t2171 * t2221;
t2216 = t2173 * t2220;
t2215 = t2175 * t2219;
t2214 = t2184 * t2230;
t2213 = t2188 * t2229;
t2212 = t2190 * t2228;
t2211 = t2192 * t2227;
t2057 = -t2084 * t2146 + t2088 * t2147;
t2061 = -t2085 * t2148 + t2089 * t2151;
t2062 = -t2086 * t2149 + t2090 * t2152;
t2063 = -t2087 * t2150 + t2091 * t2153;
t2199 = cos(qJ(1,1));
t2197 = cos(qJ(1,2));
t2195 = cos(qJ(1,3));
t2193 = sin(qJ(1,1));
t2191 = sin(qJ(1,2));
t2189 = sin(qJ(1,3));
t2187 = cos(qJ(1,4));
t2185 = sin(qJ(1,4));
t2159 = cos(t2200);
t2158 = sin(t2200);
t2093 = t2159 * g(1) + t2158 * g(2);
t2092 = t2158 * g(1) - t2159 * g(2);
t2071 = -t2087 * t2193 + t2091 * t2199;
t2070 = -t2086 * t2191 + t2090 * t2197;
t2069 = -t2085 * t2189 + t2089 * t2195;
t2068 = t2087 * t2199 + t2091 * t2193;
t2067 = t2086 * t2197 + t2090 * t2191;
t2066 = t2085 * t2195 + t2089 * t2189;
t2065 = -t2084 * t2185 + t2088 * t2187;
t2064 = t2084 * t2187 + t2088 * t2185;
t2055 = ((t2158 * t2208 + t2159 * t2204) * t2079 - 0.2e1 * (-t2158 * t2204 + t2159 * t2208) * (t2280 + (t2110 / 0.2e1 + t2109 / 0.2e1) * pkin(2))) * t2097 * t2244;
t2054 = ((t2158 * t2207 + t2159 * t2203) * t2078 - 0.2e1 * (-t2158 * t2203 + t2159 * t2207) * (t2281 + (t2108 / 0.2e1 + t2107 / 0.2e1) * pkin(2))) * t2096 * t2244;
t2053 = ((t2158 * t2206 + t2159 * t2202) * t2077 - 0.2e1 * (-t2158 * t2202 + t2159 * t2206) * (t2282 + (t2106 / 0.2e1 + t2105 / 0.2e1) * pkin(2))) * t2095 * t2244;
t2052 = ((t2158 * t2205 + t2159 * t2201) * t2073 - 0.2e1 * (-t2158 * t2201 + t2159 * t2205) * (t2283 + (t2100 / 0.2e1 + t2099 / 0.2e1) * pkin(2))) * t2094 * t2244;
t1 = [0, (t2064 * t2255 + t2066 * t2251 + t2067 * t2250 + t2068 * t2249) * t2210, (t2065 * t2255 + t2069 * t2251 + t2070 * t2250 + t2071 * t2249) * t2210, 0, (t2056 * t2255 + t2058 * t2251 + t2059 * t2250 + t2060 * t2249 + (-t2056 * t2267 - t2058 * t2263 - t2059 * t2262 - t2060 * t2261) * t2209) * t2210, (t2057 * t2255 + t2061 * t2251 + t2062 * t2250 + t2063 * t2249 + (-t2057 * t2267 - t2061 * t2263 - t2062 * t2262 - t2063 * t2261) * t2209) * t2210, 0, 0, 0, 0, 0, (t2119 * t2240 + t2125 * t2237 + t2126 * t2234 + t2127 * t2231 + (-t2073 * t2241 - t2077 * t2238 - t2078 * t2235 - t2079 * t2232) * t2209) * t2210, (-t2119 * t2222 - t2125 * t2221 - t2126 * t2220 - t2127 * t2219 + (t2073 * t2242 + t2077 * t2239 + t2078 * t2236 + t2079 * t2233) * t2209) * t2210, 0, 0, 0, -t2158 * t2092 - t2159 * t2093; 0, (t2064 * t2256 + t2066 * t2254 + t2067 * t2253 + t2068 * t2252) * t2210, (t2065 * t2256 + t2069 * t2254 + t2070 * t2253 + t2071 * t2252) * t2210, 0, (t2056 * t2256 + t2058 * t2254 + t2059 * t2253 + t2060 * t2252 + (t2056 * t2268 + t2058 * t2266 + t2059 * t2265 + t2060 * t2264) * t2209) * t2210, (t2057 * t2256 + t2061 * t2254 + t2062 * t2253 + t2063 * t2252 + (t2057 * t2268 + t2061 * t2266 + t2062 * t2265 + t2063 * t2264) * t2209) * t2210, 0, 0, 0, 0, 0, (t2118 * t2240 + t2122 * t2237 + t2123 * t2234 + t2124 * t2231 + (t2072 * t2241 + t2074 * t2238 + t2075 * t2235 + t2076 * t2232) * t2209) * t2210, (-t2118 * t2222 - t2122 * t2221 - t2123 * t2220 - t2124 * t2219 + (-t2072 * t2242 - t2074 * t2239 - t2075 * t2236 - t2076 * t2233) * t2209) * t2210, 0, 0, 0, t2159 * t2092 - t2158 * t2093; 0, (t2064 * t2226 + t2066 * t2225 + t2067 * t2224 + t2068 * t2223) * t2210, (t2065 * t2226 + t2069 * t2225 + t2070 * t2224 + t2071 * t2223) * t2210, 0, (t2218 + t2217 + t2216 + t2215 + (-t2056 * t2214 - t2058 * t2213 - t2059 * t2212 - t2060 * t2211) * t2209) * t2210, (t2057 * t2226 + t2061 * t2225 + t2062 * t2224 + t2063 * t2223 + (-t2057 * t2214 - t2061 * t2213 - t2062 * t2212 - t2063 * t2211) * t2209) * t2210, 0, 0, 0, 0, 0, (t2219 + t2220 + t2221 + t2222) * t2210 + (t2175 * (-g(3) * t2198 + t2063 * t2192) + t2173 * (-g(3) * t2196 + t2062 * t2190) + t2171 * (-g(3) * t2194 + t2061 * t2188) + t2163 * (-g(3) * t2186 + t2057 * t2184) + (-t2098 * t2218 - t2101 * t2217 - t2102 * t2216 - t2103 * t2215) * t2210) * t2209, (-t2163 * t2162 * t2292 - t2171 * t2166 * t2291 - t2173 * t2168 * t2290 - t2175 * t2170 * t2289) * t2210 + (t2175 * (g(3) * t2192 + t2063 * t2198) + t2173 * (g(3) * t2190 + t2062 * t2196) + t2171 * (g(3) * t2188 + t2061 * t2194) + t2163 * (g(3) * t2184 + t2057 * t2186) + (t2227 * t2289 + t2228 * t2290 + t2229 * t2291 + t2230 * t2292) * t2210) * t2209, 0, 0, 0, -g(3); 0, (t2064 * t2260 + t2066 * t2259 + t2067 * t2258 + t2068 * t2257) * t2210, (t2065 * t2260 + t2069 * t2259 + t2070 * t2258 + t2071 * t2257) * t2210, 0, t2052 * t2056 + t2053 * t2058 + t2054 * t2059 + t2055 * t2060 + (t2056 * t2260 + t2058 * t2259 + t2059 * t2258 + t2060 * t2257) * t2210, t2052 * t2057 + t2053 * t2061 + t2054 * t2062 + t2055 * t2063 + (t2057 * t2260 + t2061 * t2259 + t2062 * t2258 + t2063 * t2257) * t2210, 0, 0, 0, 0, 0, t2052 * t2275 + t2053 * t2273 + t2054 * t2271 + t2055 * t2269 + (t2080 * t2240 + t2081 * t2237 + t2082 * t2234 + t2083 * t2231) * t2210, -t2052 * t2276 - t2053 * t2274 - t2054 * t2272 - t2055 * t2270 + (-t2080 * t2222 - t2081 * t2221 - t2082 * t2220 - t2083 * t2219) * t2210, 0, t2092, t2093, 0;];
tau_reg  = t1;
