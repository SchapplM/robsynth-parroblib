% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G3A0
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
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:43:07
% EndTime: 2020-08-06 18:43:09
% DurationCPUTime: 2.01s
% Computational Cost: add. (912->194), mult. (1086->299), div. (111->10), fcn. (957->53), ass. (0->151)
t3180 = sin(pkin(7));
t3198 = -pkin(6) - pkin(5);
t3140 = t3198 * t3180 - pkin(1);
t3192 = cos(qJ(1,3));
t3131 = t3140 * t3192;
t3181 = cos(pkin(7));
t3186 = sin(qJ(1,3));
t3244 = t3186 * t3198;
t3245 = t3186 * t3180;
t3284 = t3131 + pkin(2) * t3245 - (pkin(2) * t3192 - t3244) * t3181;
t3194 = cos(qJ(1,2));
t3132 = t3140 * t3194;
t3188 = sin(qJ(1,2));
t3240 = t3188 * t3198;
t3241 = t3188 * t3180;
t3283 = t3132 + pkin(2) * t3241 - (pkin(2) * t3194 - t3240) * t3181;
t3196 = cos(qJ(1,1));
t3133 = t3140 * t3196;
t3190 = sin(qJ(1,1));
t3236 = t3190 * t3198;
t3237 = t3190 * t3180;
t3282 = t3133 + pkin(2) * t3237 - (pkin(2) * t3196 - t3236) * t3181;
t3269 = MDP(4) * pkin(1);
t3281 = MDP(2) / 0.2e1 + t3269 / 0.2e1;
t3279 = MDP(3) / 0.2e1;
t3278 = MDP(10) / 0.2e1;
t3277 = MDP(11) / 0.2e1;
t3272 = 0.2e1 * pkin(1);
t3271 = 0.2e1 * pkin(2);
t3270 = 0.2e1 * t3198;
t3191 = cos(qJ(3,3));
t3149 = t3191 * pkin(3) + pkin(2);
t3193 = cos(qJ(3,2));
t3150 = t3193 * pkin(3) + pkin(2);
t3195 = cos(qJ(3,1));
t3151 = t3195 * pkin(3) + pkin(2);
t3171 = qJ(1,3) + pkin(7);
t3152 = sin(t3171);
t3155 = cos(t3171);
t3185 = sin(qJ(3,3));
t3230 = pkin(7) + qJ(3,3);
t3233 = -pkin(7) + qJ(3,3);
t3268 = (t3155 * t3270 + t3186 * t3272 + t3152 * t3271 + (sin(qJ(1,3) - t3233) + sin(qJ(1,3) + t3230)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t3185 * t3271 + (sin(t3230) + sin(t3233)) * pkin(1));
t3172 = qJ(1,2) + pkin(7);
t3153 = sin(t3172);
t3156 = cos(t3172);
t3187 = sin(qJ(3,2));
t3231 = pkin(7) + qJ(3,2);
t3234 = -pkin(7) + qJ(3,2);
t3267 = (t3156 * t3270 + t3188 * t3272 + t3153 * t3271 + (sin(qJ(1,2) - t3234) + sin(qJ(1,2) + t3231)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t3187 * t3271 + (sin(t3231) + sin(t3234)) * pkin(1));
t3173 = qJ(1,1) + pkin(7);
t3154 = sin(t3173);
t3157 = cos(t3173);
t3189 = sin(qJ(3,1));
t3232 = pkin(7) + qJ(3,1);
t3235 = -pkin(7) + qJ(3,1);
t3266 = (t3157 * t3270 + t3190 * t3272 + t3154 * t3271 + (sin(qJ(1,1) - t3235) + sin(qJ(1,1) + t3232)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t3189 * t3271 + (sin(t3232) + sin(t3235)) * pkin(1));
t3182 = legFrame(3,2);
t3159 = sin(t3182);
t3162 = cos(t3182);
t3128 = t3162 * g(1) - t3159 * g(2);
t3158 = t3181 * pkin(1);
t3134 = 0.1e1 / (t3158 + t3149);
t3265 = (g(3) * t3155 + t3128 * t3152) * t3134;
t3183 = legFrame(2,2);
t3160 = sin(t3183);
t3163 = cos(t3183);
t3129 = t3163 * g(1) - t3160 * g(2);
t3135 = 0.1e1 / (t3158 + t3150);
t3264 = (g(3) * t3156 + t3129 * t3153) * t3135;
t3184 = legFrame(1,2);
t3161 = sin(t3184);
t3164 = cos(t3184);
t3130 = t3164 * g(1) - t3161 * g(2);
t3136 = 0.1e1 / (t3158 + t3151);
t3263 = (g(3) * t3157 + t3130 * t3154) * t3136;
t3142 = t3182 + t3171;
t3143 = -t3182 + t3171;
t3116 = -sin(t3142) - sin(t3143);
t3262 = t3116 * t3134;
t3144 = t3183 + t3172;
t3145 = -t3183 + t3172;
t3117 = -sin(t3144) - sin(t3145);
t3261 = t3117 * t3135;
t3146 = t3184 + t3173;
t3147 = -t3184 + t3173;
t3118 = -sin(t3146) - sin(t3147);
t3260 = t3118 * t3136;
t3119 = cos(t3143) - cos(t3142);
t3259 = t3119 * t3134;
t3120 = cos(t3145) - cos(t3144);
t3258 = t3120 * t3135;
t3121 = cos(t3147) - cos(t3146);
t3257 = t3121 * t3136;
t3256 = t3134 * t3155;
t3255 = t3134 / t3185;
t3254 = t3135 * t3156;
t3253 = t3135 / t3187;
t3252 = t3136 * t3157;
t3251 = t3136 / t3189;
t3250 = t3149 * t3192;
t3249 = t3150 * t3194;
t3248 = t3151 * t3196;
t3247 = t3185 * t3159;
t3246 = t3185 * t3162;
t3243 = t3187 * t3160;
t3242 = t3187 * t3163;
t3239 = t3189 * t3161;
t3238 = t3189 * t3164;
t3229 = pkin(3) * (-t3192 * t3181 + t3245) * t3191 ^ 2;
t3228 = pkin(3) * (-t3194 * t3181 + t3241) * t3193 ^ 2;
t3227 = pkin(3) * (-t3196 * t3181 + t3237) * t3195 ^ 2;
t3226 = ((-t3244 + t3250) * t3181 - t3131 - t3149 * t3245) * t3255;
t3225 = ((-t3240 + t3249) * t3181 - t3132 - t3150 * t3241) * t3253;
t3224 = ((-t3236 + t3248) * t3181 - t3133 - t3151 * t3237) * t3251;
t3223 = t3185 * t3265;
t3222 = t3191 * t3265;
t3221 = t3187 * t3264;
t3220 = t3193 * t3264;
t3219 = t3189 * t3263;
t3218 = t3195 * t3263;
t3125 = t3159 * g(1) + t3162 * g(2);
t3217 = t3125 * t3255;
t3126 = t3160 * g(1) + t3163 * g(2);
t3216 = t3126 * t3253;
t3127 = t3161 * g(1) + t3164 * g(2);
t3215 = t3127 * t3251;
t3214 = g(3) * t3192 + t3128 * t3186;
t3213 = g(3) * t3194 + t3129 * t3188;
t3212 = g(3) * t3196 + t3130 * t3190;
t3211 = t3159 * t3226;
t3210 = t3162 * t3226;
t3209 = t3160 * t3225;
t3208 = t3163 * t3225;
t3207 = t3161 * t3224;
t3206 = t3164 * t3224;
t3205 = -g(3) * t3152 + t3128 * t3155;
t3204 = -g(3) * t3153 + t3129 * t3156;
t3203 = -g(3) * t3154 + t3130 * t3157;
t3199 = 0.1e1 / pkin(3);
t3148 = t3158 + pkin(2);
t3115 = -g(3) * t3190 + t3130 * t3196;
t3113 = -g(3) * t3188 + t3129 * t3194;
t3111 = -g(3) * t3186 + t3128 * t3192;
t3094 = t3127 * t3189 + t3203 * t3195;
t3093 = -t3127 * t3195 + t3203 * t3189;
t3092 = t3126 * t3187 + t3204 * t3193;
t3091 = -t3126 * t3193 + t3204 * t3187;
t3090 = t3125 * t3185 + t3205 * t3191;
t3089 = -t3125 * t3191 + t3205 * t3185;
t1 = [(-(-t3164 * t3227 + (pkin(3) * t3239 - t3282 * t3164) * t3195 + t3148 * t3239) * t3215 - (-t3163 * t3228 + (pkin(3) * t3243 - t3283 * t3163) * t3193 + t3148 * t3243) * t3216 - (-t3162 * t3229 + (pkin(3) * t3247 - t3284 * t3162) * t3191 + t3148 * t3247) * t3217) * MDP(4) - g(1) * MDP(12) + ((-t3089 * t3210 - t3091 * t3208 - t3093 * t3206) * MDP(10) + (-t3090 * t3210 - t3092 * t3208 - t3094 * t3206) * MDP(11)) * t3199 + (t3111 * t3262 + t3113 * t3261 + t3115 * t3260) * t3279 + (t3116 * t3222 + t3117 * t3220 + t3118 * t3218) * t3278 + (-t3116 * t3223 - t3117 * t3221 - t3118 * t3219) * t3277 + t3281 * (t3212 * t3260 + t3213 * t3261 + t3214 * t3262); (-(t3161 * t3227 + (pkin(3) * t3238 + t3282 * t3161) * t3195 + t3148 * t3238) * t3215 - (t3160 * t3228 + (pkin(3) * t3242 + t3283 * t3160) * t3193 + t3148 * t3242) * t3216 - (t3159 * t3229 + (pkin(3) * t3246 + t3284 * t3159) * t3191 + t3148 * t3246) * t3217) * MDP(4) - g(2) * MDP(12) + ((t3089 * t3211 + t3091 * t3209 + t3093 * t3207) * MDP(10) + (t3090 * t3211 + t3092 * t3209 + t3094 * t3207) * MDP(11)) * t3199 + (t3111 * t3259 + t3113 * t3258 + t3115 * t3257) * t3279 + (t3119 * t3222 + t3120 * t3220 + t3121 * t3218) * t3278 + (-t3119 * t3223 - t3120 * t3221 - t3121 * t3219) * t3277 + t3281 * (t3212 * t3257 + t3213 * t3258 + t3214 * t3259); (-t3111 * t3256 - t3113 * t3254 - t3115 * t3252) * MDP(3) + (t3195 * ((t3151 * t3190 + t3196 * t3198) * t3181 - t3140 * t3190 + t3180 * t3248) * t3215 + t3193 * ((t3150 * t3188 + t3194 * t3198) * t3181 - t3140 * t3188 + t3180 * t3249) * t3216 + t3191 * ((t3149 * t3186 + t3192 * t3198) * t3181 - t3140 * t3186 + t3180 * t3250) * t3217) * MDP(4) + (-t3155 * t3222 - t3156 * t3220 - t3157 * t3218) * MDP(10) + (t3155 * t3223 + t3156 * t3221 + t3157 * t3219) * MDP(11) - g(3) * MDP(12) + ((t3089 * t3268 + t3091 * t3267 + t3093 * t3266) * MDP(10) + (t3090 * t3268 + t3092 * t3267 + t3094 * t3266) * MDP(11)) * t3199 + (MDP(2) + t3269) * (-t3212 * t3252 - t3213 * t3254 - t3214 * t3256);];
taugX  = t1;
