% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR6V1G2A0
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
%   see P3RPRRR6V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:37:40
% EndTime: 2020-08-06 18:37:55
% DurationCPUTime: 14.39s
% Computational Cost: add. (86064->547), mult. (101970->881), div. (13827->22), fcn. (72696->102), ass. (0->405)
t3043 = sin(pkin(7));
t3068 = -pkin(6) - pkin(5);
t2954 = t3068 * t3043;
t2940 = -t2954 + pkin(1);
t3052 = sin(qJ(1,3));
t2909 = t2940 * t3052;
t3044 = cos(pkin(7));
t3058 = cos(qJ(1,3));
t3246 = t3058 * t3068;
t3255 = t3043 * t3058;
t3358 = t2909 + pkin(2) * t3255 + (pkin(2) * t3052 + t3246) * t3044;
t3054 = sin(qJ(1,2));
t2910 = t2940 * t3054;
t3060 = cos(qJ(1,2));
t3245 = t3060 * t3068;
t3254 = t3043 * t3060;
t3357 = t2910 + pkin(2) * t3254 + (pkin(2) * t3054 + t3245) * t3044;
t3056 = sin(qJ(1,1));
t2911 = t2940 * t3056;
t3062 = cos(qJ(1,1));
t3244 = t3062 * t3068;
t3253 = t3043 * t3062;
t3356 = t2911 + pkin(2) * t3253 + (pkin(2) * t3056 + t3244) * t3044;
t3090 = pkin(2) ^ 2;
t3070 = 0.3e1 * t3090;
t3088 = pkin(3) ^ 2;
t3355 = t3070 + 0.3e1 / 0.4e1 * t3088;
t3089 = 0.1e1 / pkin(3);
t3059 = cos(qJ(3,2));
t2994 = t3059 * pkin(3);
t3065 = xDP(3);
t3207 = pkin(2) * t3044 + t2940;
t2953 = t3068 * t3044;
t3238 = pkin(2) * t3043 + t2953;
t2984 = t3044 * pkin(1);
t2951 = t2984 + pkin(2);
t3053 = sin(qJ(3,2));
t3250 = t3053 * t2951;
t3082 = 0.2e1 * qJ(3,2);
t3020 = sin(t3082);
t3306 = pkin(3) * t3020;
t3181 = t3065 * ((t2994 * t3044 + t3207) * t3060 - t3054 * (t2994 * t3043 + t3238)) / (t3306 / 0.2e1 + t3250);
t2962 = t2994 + pkin(2);
t2925 = t2984 + t2962;
t2915 = 0.1e1 / t2925;
t3036 = 0.1e1 / t3053;
t3276 = t2915 * t3036;
t3184 = ((t2962 * t3044 + t2940) * t3054 + t3060 * (t2962 * t3043 + t2953)) * t3276;
t3049 = legFrame(2,2);
t2991 = cos(t3049);
t3067 = xDP(1);
t3262 = t2991 * t3067;
t2988 = sin(t3049);
t3066 = xDP(2);
t3263 = t2988 * t3066;
t3354 = (t3181 / 0.6e1 + (t3262 / 0.6e1 - t3263 / 0.6e1) * t3184) * t3089;
t3353 = (t3181 / 0.4e1 + (t3262 / 0.4e1 - t3263 / 0.4e1) * t3184) * t3088 * t3089;
t3352 = (t3181 / 0.3e1 + (t3262 / 0.3e1 - t3263 / 0.3e1) * t3184) * t3089;
t3351 = (t3181 / 0.2e1 + (t3262 / 0.2e1 - t3263 / 0.2e1) * t3184) * t3089;
t3061 = cos(qJ(3,1));
t2995 = t3061 * pkin(3);
t2964 = t2995 + pkin(2);
t2885 = (t2964 * t3044 + t2940) * t3056 + t3062 * (t2964 * t3043 + t2953);
t3050 = legFrame(1,2);
t2989 = sin(t3050);
t2926 = t2984 + t2964;
t2918 = 0.1e1 / t2926;
t3055 = sin(qJ(3,1));
t3038 = 0.1e1 / t3055;
t3275 = t2918 * t3038;
t3183 = t2989 * t3275;
t2868 = t3066 * t2885 * t3089 * t3183;
t2992 = cos(t3050);
t3182 = t2992 * t3275;
t3248 = t3055 * t2951;
t3085 = 0.2e1 * qJ(3,1);
t3023 = sin(t3085);
t3305 = pkin(3) * t3023;
t3332 = (((t2995 * t3044 + t3207) * t3062 - t3056 * (t2995 * t3043 + t3238)) / (t3305 / 0.2e1 + t3248) * t3065 + t2885 * t3067 * t3182) * t3089;
t3350 = -0.2e1 * t2868 + 0.2e1 * t3332;
t2853 = t2868 - t3332;
t3349 = 0.2e1 * pkin(2);
t3348 = 0.4e1 * pkin(2);
t3347 = 2 * MDP(6);
t3312 = 0.2e1 * t3068;
t3346 = MDP(5) / 0.4e1;
t3345 = t2951 / 0.4e1;
t3057 = cos(qJ(3,3));
t2993 = t3057 * pkin(3);
t2961 = t2993 + pkin(2);
t2924 = t2984 + t2961;
t2912 = 0.1e1 / t2924;
t2983 = t3043 * pkin(1);
t3145 = -t3068 / 0.2e1 + t2983 / 0.2e1;
t3051 = sin(qJ(3,3));
t3048 = legFrame(3,2);
t2987 = sin(t3048);
t2990 = cos(t3048);
t3252 = t3051 * t2951;
t3034 = 0.1e1 / t3051;
t3277 = t2912 * t3034;
t3079 = 0.2e1 * qJ(3,3);
t3017 = sin(t3079);
t3307 = pkin(3) * t3017;
t2851 = (-t3065 * ((t2993 * t3044 + t3207) * t3058 - t3052 * (t2993 * t3043 + t3238)) / (t3307 / 0.2e1 + t3252) + (t2987 * t3066 - t2990 * t3067) * ((t2961 * t3044 + t2940) * t3052 + t3058 * (t2961 * t3043 + t2953)) * t3277) * t3089;
t3301 = t2851 * pkin(3);
t3213 = t3051 * t3301;
t3008 = qJ(1,3) + pkin(7);
t2944 = t3048 + t3008;
t2945 = -t3048 + t3008;
t2899 = -sin(t2944) + sin(t2945);
t2902 = cos(t2945) + cos(t2944);
t2970 = sin(t3008);
t3313 = -0.2e1 * t3065;
t2872 = t2899 * t3066 + t2902 * t3067 + t2970 * t3313;
t3338 = t2872 * t2912;
t2835 = t3145 * t3338 - t3213;
t2950 = t2983 + pkin(5);
t2943 = pkin(6) + t2950;
t2960 = t2993 + t3349;
t3077 = 0.4e1 * qJ(3,3);
t3015 = sin(t3077);
t3078 = 0.3e1 * qJ(3,3);
t3016 = sin(t3078);
t3025 = cos(t3078);
t3087 = pkin(3) * t3088;
t3076 = 0.2e1 * pkin(7);
t2996 = sin(t3076);
t2997 = cos(t3076);
t3091 = pkin(1) ^ 2;
t3042 = t3068 ^ 2;
t3237 = t3042 / 0.2e1 + t3091;
t3148 = 0.3e1 / 0.8e1 * t3088 + t3090 / 0.2e1 + t3237;
t3215 = pkin(1) * t2954;
t3235 = t3042 + t3091;
t3111 = -0.8e1 * (t3235 + t3355) * t2984 - 0.8e1 * (t3148 - t3215) * t3349 + 0.8e1 * (-pkin(2) * t2997 + t2996 * t3068) * t3091;
t3116 = pkin(3) * (-pkin(2) * t3057 + t2951 * t3025);
t3214 = pkin(1) * t2953;
t3239 = t3091 * t2996;
t3218 = 0.2e1 * t3239;
t3120 = -0.4e1 * t3214 + t3218;
t2965 = t3091 * t2997;
t3331 = -t2965 - t3088 / 0.2e1;
t3178 = -0.2e1 * t3090 + t3331;
t3147 = -t3091 + t3178;
t3075 = -0.2e1 * t3091;
t3149 = -0.2e1 * t2965 - 0.4e1 * t3090 + t3075 - t3088;
t3169 = -0.2e1 * t3214 + t3239;
t3180 = -t3089 / 0.4e1;
t3285 = t2872 / 0.2e1;
t3190 = t2912 * t3285;
t3296 = 0.2e1 * pkin(3);
t3203 = t2943 * t3296;
t3216 = t2960 * t2984;
t3311 = -0.6e1 * t3088;
t3219 = t2951 * t3311;
t3223 = -0.2e1 * t2943 * t3088;
t2941 = -0.2e1 * t3215;
t3074 = 0.2e1 * t3091;
t3224 = pkin(2) * t2984;
t3230 = -0.4e1 * pkin(3) * (0.6e1 * t3224 + t2941 + t3074 + t3070 + t3042 - t3331);
t3243 = t3088 * cos(t3077);
t3071 = 0.2e1 * t3090;
t3234 = t3071 + t3091;
t2905 = 0.4e1 * t3224 + t2965 + t3234;
t3026 = cos(t3079);
t3278 = t2905 * t3026;
t3308 = pkin(2) * t2943;
t2894 = t3169 + 0.2e1 * t3308;
t3279 = t2894 * t3026;
t3295 = 0.4e1 * pkin(3);
t3315 = -0.2e1 * t3088 - 0.2e1 * t2905;
t2822 = ((t3025 * t3223 + (t2943 * t2960 + t3169 - t3279) * t3296) * t2851 + (-t3015 * t3087 + t3016 * t3219 + t3017 * t3230 + t3051 * t3111) * t3190) / (t3278 + t3243 / 0.2e1 - 0.2e1 * t3216 + 0.2e1 * t3116 + t3147) * t3180 * t3338 - (t2835 * t3348 + (t3017 * t3315 - t3088 * t3015 + (-t2951 * t3016 - t2984 * t3051) * t3295) * t2851 + (-0.2e1 * t3279 + (-t3025 + t3057) * t3203 + t3120) * t3190) / (0.4e1 * t3116 + t3149 - 0.4e1 * t3216 + t3243 + 0.2e1 * t3278) * t2851;
t2913 = 0.1e1 / t2924 ^ 2;
t2832 = (-t2835 + t3213) * t2913 * t3285;
t2848 = t2851 ^ 2;
t3158 = t3051 * t3190;
t3176 = MDP(4) * t3091 + MDP(1);
t3189 = t3057 * t3285;
t3222 = 0.2e1 * t2951 * t3057;
t3251 = t3051 * t3057;
t3271 = t2950 * t3057;
t3274 = t2950 * t3051;
t3291 = t2848 * t3057;
t3292 = t2822 * t3057;
t3318 = 0.2e1 * t2851;
t3326 = MDP(10) * (t2822 * t3274 - t2832 * t3222 + (t2851 * t3271 / 0.2e1 + t2951 * t3158) * t3318) + MDP(11) * (t2822 * t3271 + 0.2e1 * t2832 * t3252 + (t2912 * t2951 * t3189 - t2851 * t3274 / 0.2e1) * t3318) - MDP(7) * (t2822 * t3051 + t3291) - MDP(8) * (-t2848 * t3051 + t3292) - t2832 * (MDP(5) * t3051 ^ 2 + t3251 * t3347 + t3176);
t3343 = -t3326 * t2912 / 0.2e1;
t3299 = t2853 * pkin(3);
t3211 = t3055 * t3299;
t3014 = qJ(1,1) + pkin(7);
t2948 = t3050 + t3014;
t2949 = -t3050 + t3014;
t2901 = -sin(t2948) + sin(t2949);
t2904 = cos(t2949) + cos(t2948);
t2972 = sin(t3014);
t2874 = t2901 * t3066 + t2904 * t3067 + t2972 * t3313;
t3336 = t2874 * t2918;
t2837 = t3145 * t3336 - t3211;
t3283 = t2874 / 0.2e1;
t2867 = t2918 * t3283;
t2846 = t2867 + t3350;
t2847 = t2867 - t3350;
t3002 = pkin(7) + qJ(3,1);
t2956 = 0.2e1 * t3002;
t3007 = -pkin(7) + qJ(3,1);
t2958 = 0.2e1 * t3007;
t2963 = t2995 + t3349;
t2967 = sin(t3002);
t2969 = sin(t3007);
t3001 = pkin(7) + t3085;
t3006 = -pkin(7) + t3085;
t3084 = 0.3e1 * qJ(3,1);
t3012 = t3084 + pkin(7);
t3013 = t3084 - pkin(7);
t3083 = 0.4e1 * qJ(3,1);
t3021 = sin(t3083);
t3022 = sin(t3084);
t3031 = cos(t3084);
t3032 = cos(t3085);
t3177 = pkin(3) * (t3031 - t3061);
t3314 = -0.2e1 * t3044;
t3179 = cos(t3001) + cos(t3006) + t3314;
t3241 = t3088 * cos(t3083);
t3310 = pkin(2) * pkin(3);
t2819 = ((t3031 * t3223 + (-t2894 * t3032 + t2943 * t2963 + t3169) * t3296) * t2853 + (-t3021 * t3087 + t3022 * t3219 + t3023 * t3230 + t3055 * t3111) * t2867) / (t2905 * t3032 + t3241 / 0.2e1 - 0.2e1 * t2963 * t2984 + (-t3061 * pkin(2) + t2951 * t3031) * t3296 + t3147) * t3180 * t3336 - (t2837 * t3348 + ((t2867 - t2853) * sin(t2958) - (t2867 + t2853) * sin(t2956)) * t3091 + (-0.2e1 * (t3088 + t3234) * t3023 - 0.4e1 * t3022 * t3310 - t3088 * t3021) * t2853 + (t3218 + (t3032 * t3348 + 0.2e1 * t3177) * t3068) * t2867 + ((sin(t3006) * t2846 - sin(t3001) * t2847) * t3349 + t3179 * t2867 * t3312 + ((-sin(t3012) - t2969) * t2847 + (sin(t3013) + t2967) * t2846) * pkin(3)) * pkin(1)) / ((t3074 + 0.4e1 * t3090) * t3032 + t3241 + (cos(t2958) + cos(t2956)) * t3091 + t3177 * t3348 + (t3179 * t3348 + (cos(t3013) + cos(t3012) - cos(t3007) - cos(t3002)) * t3296) * pkin(1) + t3149) * t2853;
t2919 = 0.1e1 / t2926 ^ 2;
t3187 = t2919 * t3283;
t2834 = (-t2837 + t3211) * t3187;
t2850 = t2853 ^ 2;
t3156 = t3055 * t2867;
t3267 = t2951 * t3061;
t3220 = 0.2e1 * t3267;
t3247 = t3055 * t3061;
t3269 = t2950 * t3061;
t3272 = t2950 * t3055;
t3289 = t2850 * t3061;
t3293 = t2819 * t3061;
t3316 = 0.2e1 * t2853;
t3324 = MDP(10) * (t2819 * t3272 - t2834 * t3220 + (t2853 * t3269 / 0.2e1 + t2951 * t3156) * t3316) + MDP(11) * (t2819 * t3269 + 0.2e1 * t2834 * t3248 + (t2867 * t3267 - t2853 * t3272 / 0.2e1) * t3316) - MDP(7) * (t2819 * t3055 + t3289) - MDP(8) * (-t2850 * t3055 + t3293) - t2834 * (MDP(5) * t3055 ^ 2 + t3247 * t3347 + t3176);
t3342 = -t3324 * t2918 / 0.2e1;
t2852 = (-t3181 + (-t3262 + t3263) * t3184) * t3089;
t3300 = t2852 * pkin(3);
t3212 = t3053 * t3300;
t3011 = qJ(1,2) + pkin(7);
t2946 = t3049 + t3011;
t2947 = -t3049 + t3011;
t2900 = -sin(t2946) + sin(t2947);
t2903 = cos(t2947) + cos(t2946);
t2971 = sin(t3011);
t2873 = t2900 * t3066 + t2903 * t3067 + t2971 * t3313;
t3337 = t2873 * t2915;
t2836 = t3145 * t3337 - t3212;
t3284 = t2873 / 0.2e1;
t2866 = t2915 * t3284;
t2865 = t3091 * t2866;
t3000 = pkin(7) + qJ(3,2);
t2955 = 0.2e1 * t3000;
t3005 = -pkin(7) + qJ(3,2);
t2957 = 0.2e1 * t3005;
t2966 = sin(t3000);
t2968 = sin(t3005);
t2999 = pkin(7) + t3082;
t2973 = cos(t2999);
t2974 = cos(t3000);
t3004 = -pkin(7) + t3082;
t2976 = cos(t3004);
t2977 = cos(t3005);
t3081 = 0.3e1 * qJ(3,2);
t2998 = pkin(7) + t3081;
t3003 = -pkin(7) + t3081;
t3009 = qJ(3,2) + t3076;
t3010 = qJ(3,2) - 0.2e1 * pkin(7);
t3080 = 0.4e1 * qJ(3,2);
t3018 = sin(t3080);
t3019 = sin(t3081);
t3028 = cos(t3081);
t3029 = cos(t3082);
t3146 = t2865 + (t3042 + t3355) * t2866;
t3157 = t3053 * t2866;
t3309 = pkin(1) * t3068;
t3225 = pkin(3) * t3309;
t3228 = 0.4e1 * t3300;
t3229 = pkin(3) * t3075;
t3236 = t3028 - t3059;
t3242 = t3088 * cos(t3080);
t2816 = (t3088 * t2852 * t3028 * t3312 + t3218 * t3300 + (t3308 + (-t2984 - t2994 / 0.2e1) * t3068) * t3228 + (-0.8e1 * (t3088 / 0.4e1 + 0.3e1 / 0.2e1 * t3090 + t3237) * t3306 - t3087 * t3018) * t2866 + 0.4e1 * (cos(t3010) - cos(t3009)) * t3068 * t2865 + (sin(t2955) * t3229 + 0.4e1 * t2976 * t3225) * (t2866 - t3351) + (sin(t2957) * t3229 - 0.4e1 * t2973 * t3225) * (t2866 + t3351) + (-0.4e1 * (t3146 - t3353) * t2968 - 0.4e1 * (t3146 + t3353) * t2966 - 0.3e1 * ((t2866 + t3352) * sin(t3003) + (t2866 - t3352) * sin(t2998)) * t3088 - 0.12e2 * ((t2866 + t3354) * sin(t3004) + (t2866 - t3354) * sin(t2999)) * t3310) * pkin(1) + ((t3019 * t3311 + 0.8e1 * (t2977 - t2974) * t3309) * t2866 - 0.4e1 * (sin(t3010) + sin(t3009)) * t2865 + t3029 * t3068 * t3228 - 0.16e2 * t3148 * t3157) * pkin(2)) / (t3071 * t3029 + t3242 / 0.2e1 + (cos(t2957) / 0.2e1 + cos(t2955) / 0.2e1 + t3029 - 0.1e1) * t3091 + t3236 * pkin(2) * t3296 + ((t2973 + t2976 + t3314) * t3349 + (cos(t3003) + cos(t2998) - t2974 - t2977) * pkin(3)) * pkin(1) + t3178) * t3180 * t3337 - (t2836 * t3348 + (t3020 * t3315 - t3088 * t3018 + (-t2951 * t3019 - t2984 * t3053) * t3295) * t2852 + (-0.2e1 * t2894 * t3029 - t3203 * t3236 + t3120) * t2866) / (0.2e1 * t2905 * t3029 + t3242 - 0.4e1 * (t2994 + t3349) * t2984 + (-t3059 * pkin(2) + t2951 * t3028) * t3295 + t3149) * t2852;
t2916 = 0.1e1 / t2925 ^ 2;
t3188 = t2916 * t3284;
t2833 = (-t2836 + t3212) * t3188;
t2849 = t2852 ^ 2;
t3268 = t2951 * t3059;
t3221 = 0.2e1 * t3268;
t3249 = t3053 * t3059;
t3270 = t2950 * t3059;
t3273 = t2950 * t3053;
t3290 = t2849 * t3059;
t3294 = t2816 * t3059;
t3317 = 0.2e1 * t2852;
t3325 = MDP(10) * (t2816 * t3273 - t2833 * t3221 + (t2852 * t3270 / 0.2e1 + t2951 * t3157) * t3317) + MDP(11) * (t2816 * t3270 + 0.2e1 * t2833 * t3250 + (t2866 * t3268 - t2852 * t3273 / 0.2e1) * t3317) - MDP(7) * (t2816 * t3053 + t3290) - MDP(8) * (-t2849 * t3053 + t3294) - t2833 * (MDP(5) * t3053 ^ 2 + t3249 * t3347 + t3176);
t3341 = -t3325 * t2915 / 0.2e1;
t3340 = t3347 / 0.4e1;
t3227 = 0.2e1 * MDP(5);
t3339 = t3227 / 0.4e1;
t3288 = t2872 ^ 2 * t2913;
t3287 = t2873 ^ 2 * t2916;
t3286 = t2874 ^ 2 * t2919;
t3266 = t2961 * t3052;
t2877 = (t3246 + t3266) * t3044 + t2909 + t2961 * t3255;
t3335 = t2877 * t2912;
t3265 = t2962 * t3054;
t2879 = (t3245 + t3265) * t3044 + t2910 + t2962 * t3254;
t3334 = t2879 * t2915;
t3264 = t2964 * t3056;
t2881 = (t3244 + t3264) * t3044 + t2911 + t2964 * t3253;
t3333 = t2881 * t2918;
t3322 = -0.2e1 * pkin(1);
t3321 = -0.2e1 * pkin(2);
t3039 = t3057 ^ 2;
t2979 = 0.2e1 * t3039 - 0.1e1;
t3040 = t3059 ^ 2;
t2980 = 0.2e1 * t3040 - 0.1e1;
t3041 = t3061 ^ 2;
t2981 = 0.2e1 * t3041 - 0.1e1;
t3304 = pkin(3) * t3051;
t3303 = pkin(3) * t3053;
t3302 = pkin(3) * t3055;
t3232 = pkin(7) + qJ(3,3);
t3233 = -pkin(7) + qJ(3,3);
t3282 = (t2970 * t3312 + cos(t3008) * t3321 + t3058 * t3322 + (-cos(qJ(1,3) - t3233) - cos(qJ(1,3) + t3232)) * pkin(3)) / (t3051 * t3349 + t3307 + (sin(t3232) + sin(t3233)) * pkin(1));
t3281 = (t2971 * t3312 + cos(t3011) * t3321 + t3060 * t3322 + (-cos(qJ(1,2) - t3005) - cos(qJ(1,2) + t3000)) * pkin(3)) / (t3053 * t3349 + t3306 + (t2966 + t2968) * pkin(1));
t3280 = (t2972 * t3312 + cos(t3014) * t3321 + t3062 * t3322 + (-cos(qJ(1,1) - t3007) - cos(qJ(1,1) + t3002)) * pkin(3)) / (t3055 * t3349 + t3305 + (t2967 + t2969) * pkin(1));
t3261 = t3034 * t3039;
t3260 = t3034 * t3057;
t3259 = t3036 * t3040;
t3258 = t3036 * t3059;
t3257 = t3038 * t3041;
t3256 = t3038 * t3061;
t3210 = pkin(3) * (t3044 * t3052 + t3255) * t3039;
t3209 = pkin(3) * (t3044 * t3056 + t3253) * t3041;
t3208 = (t3044 * t3054 + t3254) * t3040 * pkin(3);
t2895 = t2941 + t3090 + 0.2e1 * t3224 + t3235;
t2829 = ((-t2943 * t3213 + (pkin(3) * t3222 + t3039 * t3088 + t2895) * t3190) * t2913 * t3189 - ((t2943 * t3158 - t3301) * t3057 - t2851 * t2951) * t2912 * t3301) * t3034;
t3206 = t2829 * t3034 * MDP(4);
t2830 = ((-t2943 * t3212 + (pkin(3) * t3221 + t3040 * t3088 + t2895) * t2866) * t3059 * t3188 - ((t2943 * t3157 - t3300) * t3059 - t2852 * t2951) * t2915 * t3300) * t3036;
t3205 = t2830 * t3036 * MDP(4);
t2831 = ((-t2943 * t3211 + (pkin(3) * t3220 + t3041 * t3088 + t2895) * t2867) * t3061 * t3187 - ((t2943 * t3156 - t3299) * t3061 - t2853 * t2951) * t2918 * t3299) * t3038;
t3204 = t2831 * t3038 * MDP(4);
t3202 = t2832 * t3335;
t3201 = t2832 * t3282;
t3200 = t2833 * t3334;
t3199 = t2833 * t3281;
t3198 = t2834 * t3333;
t3197 = t2834 * t3280;
t3196 = t2851 * t2872 * t2913;
t3195 = t2852 * t2873 * t2916;
t3194 = t2853 * t2874 * t2919;
t3193 = t3335 * t3288;
t3192 = t3334 * t3287;
t3191 = t3333 * t3286;
t3186 = t2877 * t3277;
t3185 = t2879 * t3276;
t3167 = t2979 * t3196;
t3166 = t2980 * t3195;
t3165 = t2981 * t3194;
t3164 = t3057 * t3193;
t3163 = t3282 * t3288;
t3162 = t3059 * t3192;
t3161 = t3281 * t3287;
t3160 = t3061 * t3191;
t3159 = t3280 * t3286;
t3155 = t2987 * t3186;
t3154 = t2990 * t3186;
t3153 = t2988 * t3185;
t3152 = t2991 * t3185;
t3151 = t2881 * t3183;
t3150 = t2881 * t3182;
t3144 = t2832 * t3057 * t3186;
t3143 = t2833 * t3059 * t3185;
t3142 = t3198 * t3256;
t3141 = t3196 * t3251;
t3140 = t3195 * t3249;
t3139 = t3194 * t3247;
t3138 = t3057 * t3163;
t3137 = t3059 * t3161;
t3136 = t3061 * t3159;
t3134 = t2987 * t3164;
t3133 = t2988 * t3162;
t3132 = t2989 * t3160;
t3131 = t2990 * t3164;
t3130 = t2991 * t3162;
t3129 = t2992 * t3160;
t3104 = t2912 * (t3206 + (t2822 * t3260 - t2848) * MDP(10) + (-t2848 * t3260 - t2822) * MDP(11));
t3103 = t2918 * (t3204 + (t2819 * t3256 - t2850) * MDP(10) + (-t2850 * t3256 - t2819) * MDP(11));
t3102 = (t3205 + (t2816 * t3258 - t2849) * MDP(10) + (-t2849 * t3258 - t2816) * MDP(11)) * t2915;
t2856 = -t2981 * t3286 / 0.4e1;
t2855 = -t2980 * t3287 / 0.4e1;
t2854 = -t2979 * t3288 / 0.4e1;
t2828 = -t2831 * t3055 - t2834 * t3269;
t2827 = t2831 * t3061 - t2834 * t3272;
t2826 = -t2830 * t3053 - t2833 * t3270;
t2825 = t2830 * t3059 - t2833 * t3273;
t2824 = -t2829 * t3051 - t2832 * t3271;
t2823 = t2829 * t3057 - t2832 * t3274;
t1 = [(t2992 * t3209 + (t2989 * t3302 + t2992 * t3356) * t3061 + t2989 * t3248) * t3103 + (t2991 * t3208 + (t2988 * t3303 + t2991 * t3357) * t3059 + t2988 * t3250) * t3102 + (t2990 * t3210 + (t2987 * t3304 + t2990 * t3358) * t3057 + t2987 * t3252) * t3104 + (t2902 * t3141 + t2903 * t3140 + t2904 * t3139) * t3339 + (t2902 * t3167 + t2903 * t3166 + t2904 * t3165) * t3340 + ((-t2854 * t3154 - t2855 * t3152 - t2856 * t3150) * MDP(6) + (-t2990 * t3202 - t2991 * t3200 - t2992 * t3198) * MDP(7) + (-t2990 * t3144 - t2991 * t3143 - t2992 * t3142) * MDP(8) + (-t2816 * t3152 - t2819 * t3150 - t2822 * t3154) * MDP(9) + (-t2823 * t3154 - t2825 * t3152 - t2827 * t3150) * MDP(10) + (-t2824 * t3154 - t2826 * t3152 - t2828 * t3150) * MDP(11) + (t3129 + t3130 + t3131) * t3346 + ((-t2990 * t3193 - t2991 * t3192 - t2992 * t3191) * MDP(10) + (-t3034 * t3131 - t3036 * t3130 - t3038 * t3129) * MDP(11)) * t3345) * t3089 + t2902 * t3343 + t2903 * t3341 + t2904 * t3342; (-t2989 * t3209 + (-t2989 * t3356 + t2992 * t3302) * t3061 + t2992 * t3248) * t3103 + (-t2988 * t3208 + (-t2988 * t3357 + t2991 * t3303) * t3059 + t2991 * t3250) * t3102 + (-t2987 * t3210 + (-t2987 * t3358 + t2990 * t3304) * t3057 + t2990 * t3252) * t3104 + (t2899 * t3141 + t2900 * t3140 + t2901 * t3139) * t3339 + (t2899 * t3167 + t2900 * t3166 + t2901 * t3165) * t3340 + ((t2854 * t3155 + t2855 * t3153 + t2856 * t3151) * MDP(6) + (t2987 * t3202 + t2988 * t3200 + t2989 * t3198) * MDP(7) + (t2987 * t3144 + t2988 * t3143 + t2989 * t3142) * MDP(8) + (t2816 * t3153 + t2819 * t3151 + t2822 * t3155) * MDP(9) + (t2823 * t3155 + t2825 * t3153 + t2827 * t3151) * MDP(10) + (t2824 * t3155 + t2826 * t3153 + t2828 * t3151) * MDP(11) + (-t3132 - t3133 - t3134) * t3346 + ((t2987 * t3193 + t2988 * t3192 + t2989 * t3191) * MDP(10) + (t3034 * t3134 + t3036 * t3133 + t3038 * t3132) * MDP(11)) * t3345) * t3089 + t2899 * t3343 + t2900 * t3341 + t2901 * t3342; (-t2970 * t3141 - t2971 * t3140 - t2972 * t3139) * t3227 / 0.2e1 + (-t2970 * t3167 - t2971 * t3166 - t2972 * t3165) * t3347 / 0.2e1 + ((t3061 * t3204 + (t2819 * t3257 - t3289) * MDP(10) + (-t2850 * t3257 - t3293) * MDP(11)) * ((t2964 * t3062 - t3056 * t3068) * t3044 + t2940 * t3062 - t3043 * t3264) + t3324 * t2972) * t2918 + ((t3059 * t3205 + (t2816 * t3259 - t3290) * MDP(10) + (-t2849 * t3259 - t3294) * MDP(11)) * ((t2962 * t3060 - t3054 * t3068) * t3044 + t2940 * t3060 - t3043 * t3265) + t3325 * t2971) * t2915 + ((t3057 * t3206 + (t2822 * t3261 - t3291) * MDP(10) + (-t2848 * t3261 - t3292) * MDP(11)) * ((t2961 * t3058 - t3052 * t3068) * t3044 + t2940 * t3058 - t3043 * t3266) + t3326 * t2970) * t2912 + ((t2854 * t3282 + t2855 * t3281 + t2856 * t3280) * MDP(6) + (t3051 * t3201 + t3053 * t3199 + t3055 * t3197) * MDP(7) + (t3057 * t3201 + t3059 * t3199 + t3061 * t3197) * MDP(8) + (t2816 * t3281 + t2819 * t3280 + t2822 * t3282) * MDP(9) + (t2823 * t3282 + t2825 * t3281 + t2827 * t3280) * MDP(10) + (t2824 * t3282 + t2826 * t3281 + t2828 * t3280) * MDP(11) + (-t3051 * t3138 - t3053 * t3137 - t3055 * t3136) * t3346 + ((t3051 * t3163 + t3053 * t3161 + t3055 * t3159) * MDP(10) + (t3136 + t3137 + t3138) * MDP(11)) * t3345) * t3089;];
taucX  = t1;
