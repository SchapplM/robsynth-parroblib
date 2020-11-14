% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G1A0
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
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:48:29
% EndTime: 2020-08-06 18:48:31
% DurationCPUTime: 1.49s
% Computational Cost: add. (750->146), mult. (942->259), div. (117->7), fcn. (954->29), ass. (0->99)
t3305 = MDP(3) - MDP(6);
t3255 = cos(pkin(7));
t3301 = -MDP(2) - MDP(4) * t3255 + MDP(5) * sin(pkin(7));
t3300 = 2 * pkin(1);
t3299 = -pkin(5) - pkin(6);
t3251 = pkin(7) + qJ(3,3);
t3236 = cos(t3251);
t3298 = pkin(3) * t3236;
t3252 = pkin(7) + qJ(3,2);
t3237 = cos(t3252);
t3297 = pkin(3) * t3237;
t3253 = pkin(7) + qJ(3,1);
t3238 = cos(t3253);
t3296 = pkin(3) * t3238;
t3248 = -qJ(2,3) + t3299;
t3245 = 0.1e1 / t3248;
t3256 = legFrame(3,3);
t3239 = sin(t3256);
t3242 = cos(t3256);
t3259 = sin(qJ(1,3));
t3262 = cos(qJ(1,3));
t3272 = t3239 * t3259 - t3242 * t3262;
t3293 = t3272 * t3245;
t3249 = -qJ(2,2) + t3299;
t3246 = 0.1e1 / t3249;
t3257 = legFrame(2,3);
t3240 = sin(t3257);
t3243 = cos(t3257);
t3260 = sin(qJ(1,2));
t3263 = cos(qJ(1,2));
t3271 = t3240 * t3260 - t3243 * t3263;
t3292 = t3271 * t3246;
t3250 = -qJ(2,1) + t3299;
t3247 = 0.1e1 / t3250;
t3258 = legFrame(1,3);
t3241 = sin(t3258);
t3244 = cos(t3258);
t3261 = sin(qJ(1,1));
t3264 = cos(qJ(1,1));
t3270 = t3241 * t3261 - t3244 * t3264;
t3291 = t3270 * t3247;
t3220 = t3239 * t3262 + t3242 * t3259;
t3290 = t3220 * t3245;
t3221 = t3240 * t3263 + t3243 * t3260;
t3289 = t3221 * t3246;
t3222 = t3241 * t3264 + t3244 * t3261;
t3288 = t3222 * t3247;
t3193 = t3220 * g(1) + t3272 * g(2);
t3287 = t3245 * t3193;
t3233 = sin(t3251);
t3286 = t3245 * t3233;
t3194 = t3221 * g(1) + t3271 * g(2);
t3285 = t3246 * t3194;
t3234 = sin(t3252);
t3284 = t3246 * t3234;
t3195 = t3222 * g(1) + t3270 * g(2);
t3283 = t3247 * t3195;
t3235 = sin(t3253);
t3282 = t3247 * t3235;
t3281 = t3236 * t3287;
t3280 = t3237 * t3285;
t3279 = t3238 * t3283;
t3278 = t3193 * t3286;
t3277 = t3194 * t3284;
t3276 = t3195 * t3282;
t3230 = 0.1e1 / t3236;
t3275 = t3230 * t3286;
t3231 = 0.1e1 / t3237;
t3274 = t3231 * t3284;
t3232 = 0.1e1 / t3238;
t3273 = t3232 * t3282;
t3223 = t3239 * g(1) - t3242 * g(2);
t3226 = t3242 * g(1) + t3239 * g(2);
t3197 = t3223 * t3262 + t3226 * t3259;
t3196 = t3223 * t3259 - t3226 * t3262;
t3224 = t3240 * g(1) - t3243 * g(2);
t3227 = t3243 * g(1) + t3240 * g(2);
t3199 = t3224 * t3263 + t3227 * t3260;
t3198 = t3224 * t3260 - t3227 * t3263;
t3225 = t3241 * g(1) - t3244 * g(2);
t3228 = t3244 * g(1) + t3241 * g(2);
t3201 = t3225 * t3264 + t3228 * t3261;
t3200 = t3225 * t3261 - t3228 * t3264;
t3266 = 0.1e1 / pkin(3);
t3265 = 0.2e1 * pkin(7);
t3229 = t3255 * pkin(2) + pkin(1);
t3216 = t3229 * t3264 - t3261 * t3250;
t3215 = t3229 * t3263 - t3260 * t3249;
t3214 = t3229 * t3262 - t3259 * t3248;
t3213 = t3261 * t3229 + t3264 * t3250;
t3212 = t3260 * t3229 + t3263 * t3249;
t3211 = t3259 * t3229 + t3262 * t3248;
t3210 = (g(1) * t3264 + g(2) * t3261) * t3244 - (g(1) * t3261 - g(2) * t3264) * t3241;
t3209 = (t3263 * g(1) + t3260 * g(2)) * t3243 - (t3260 * g(1) - t3263 * g(2)) * t3240;
t3208 = (t3262 * g(1) + t3259 * g(2)) * t3242 - (t3259 * g(1) - t3262 * g(2)) * t3239;
t3192 = -t3228 * (-t3261 * pkin(1) + t3264 * qJ(2,1)) + t3225 * (t3264 * pkin(1) + t3261 * qJ(2,1));
t3191 = -t3227 * (-t3260 * pkin(1) + t3263 * qJ(2,2)) + t3224 * (t3263 * pkin(1) + t3260 * qJ(2,2));
t3190 = -t3226 * (-t3259 * pkin(1) + t3262 * qJ(2,3)) + t3223 * (t3262 * pkin(1) + t3259 * qJ(2,3));
t1 = [(-(-t3270 * t3192 - (-t3213 * t3241 + t3216 * t3244 - t3270 * t3296) * t3201) * t3247 - (-t3271 * t3191 - (-t3212 * t3240 + t3215 * t3243 - t3271 * t3297) * t3199) * t3246 - (-t3272 * t3190 - (-t3211 * t3239 + t3214 * t3242 - t3272 * t3298) * t3197) * t3245) * MDP(7) + (t3270 * t3279 + t3271 * t3280 + t3272 * t3281) * MDP(13) + (-t3270 * t3276 - t3271 * t3277 - t3272 * t3278) * MDP(14) - g(1) * MDP(15) + t3301 * (-t3197 * t3293 - t3199 * t3292 - t3201 * t3291) - t3305 * (t3196 * t3293 + t3198 * t3292 + t3200 * t3291); (-(t3222 * t3192 - (t3213 * t3244 + t3216 * t3241 + t3222 * t3296) * t3201) * t3247 - (t3221 * t3191 - (t3212 * t3243 + t3215 * t3240 + t3221 * t3297) * t3199) * t3246 - (t3220 * t3190 - (t3211 * t3242 + t3214 * t3239 + t3220 * t3298) * t3197) * t3245) * MDP(7) + (-t3220 * t3281 - t3221 * t3280 - t3222 * t3279) * MDP(13) + (t3220 * t3278 + t3221 * t3277 + t3222 * t3276) * MDP(14) - g(2) * MDP(15) + t3301 * (t3197 * t3290 + t3199 * t3289 + t3201 * t3288) + t3305 * (t3196 * t3290 + t3198 * t3289 + t3200 * t3288); (-t3192 * t3273 + t3247 * (pkin(3) * sin(0.2e1 * t3253) + t3235 * t3300 + (sin(t3265 + qJ(3,1)) + sin(qJ(3,1))) * pkin(2)) * t3232 * t3201 / 0.2e1 - t3191 * t3274 + t3246 * (pkin(3) * sin(0.2e1 * t3252) + t3234 * t3300 + (sin(qJ(3,2)) + sin(t3265 + qJ(3,2))) * pkin(2)) * t3231 * t3199 / 0.2e1 - t3190 * t3275 + t3245 * (pkin(3) * sin(0.2e1 * t3251) + t3233 * t3300 + (sin(t3265 + qJ(3,3)) + sin(qJ(3,3))) * pkin(2)) * t3230 * t3197 / 0.2e1) * MDP(7) + (-t3278 - t3277 - t3276 + (t3232 * (-g(3) * t3238 + t3210 * t3235) + t3231 * (-g(3) * t3237 + t3209 * t3234) + t3230 * (-g(3) * t3236 + t3208 * t3233)) * t3266) * MDP(13) + (t3233 ^ 2 * t3230 * t3287 + t3234 ^ 2 * t3231 * t3285 + t3235 ^ 2 * t3232 * t3283 + (t3232 * (g(3) * t3235 + t3210 * t3238) + t3231 * (g(3) * t3234 + t3209 * t3237) + t3230 * (g(3) * t3233 + t3208 * t3236)) * t3266) * MDP(14) - g(3) * MDP(15) + t3301 * (t3197 * t3275 + t3199 * t3274 + t3201 * t3273) + t3305 * (t3196 * t3275 + t3198 * t3274 + t3200 * t3273);];
taugX  = t1;
