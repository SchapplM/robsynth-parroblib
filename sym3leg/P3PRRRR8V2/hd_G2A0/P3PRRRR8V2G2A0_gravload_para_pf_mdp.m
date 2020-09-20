% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G2A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:50:52
% EndTime: 2020-08-06 17:50:55
% DurationCPUTime: 2.29s
% Computational Cost: add. (1098->222), mult. (2547->446), div. (81->4), fcn. (2667->22), ass. (0->179)
t3100 = sin(qJ(2,3));
t3106 = cos(qJ(2,3));
t3111 = pkin(7) + pkin(6);
t3074 = pkin(2) * t3100 - t3111 * t3106;
t3093 = sin(pkin(4));
t3095 = cos(pkin(4));
t3099 = sin(qJ(3,3));
t3172 = t3095 * t3099;
t3045 = pkin(3) * t3172 + t3074 * t3093;
t3105 = cos(qJ(3,3));
t3176 = t3093 * t3100;
t3199 = pkin(3) * t3105 ^ 2;
t3030 = 0.1e1 / (pkin(2) * t3172 + t3045 * t3105 + t3176 * t3199);
t3094 = cos(pkin(8));
t3173 = t3094 * t3095;
t3066 = -t3093 * g(1) + g(2) * t3173;
t3067 = g(1) * t3173 + t3093 * g(2);
t3096 = legFrame(3,2);
t3083 = sin(t3096);
t3086 = cos(t3096);
t3071 = t3086 * g(1) - t3083 * g(2);
t3092 = sin(pkin(8));
t3178 = t3092 * t3095;
t3080 = g(3) * t3178;
t3082 = g(3) * t3094;
t3208 = t3030 * ((t3066 * t3083 - t3067 * t3086 + t3080) * t3106 + t3100 * (t3071 * t3092 + t3082));
t3102 = sin(qJ(2,2));
t3108 = cos(qJ(2,2));
t3075 = pkin(2) * t3102 - t3111 * t3108;
t3101 = sin(qJ(3,2));
t3170 = t3095 * t3101;
t3046 = pkin(3) * t3170 + t3075 * t3093;
t3107 = cos(qJ(3,2));
t3175 = t3093 * t3102;
t3198 = pkin(3) * t3107 ^ 2;
t3031 = 0.1e1 / (pkin(2) * t3170 + t3046 * t3107 + t3175 * t3198);
t3097 = legFrame(2,2);
t3084 = sin(t3097);
t3087 = cos(t3097);
t3072 = t3087 * g(1) - t3084 * g(2);
t3207 = t3031 * ((t3066 * t3084 - t3067 * t3087 + t3080) * t3108 + t3102 * (t3072 * t3092 + t3082));
t3104 = sin(qJ(2,1));
t3110 = cos(qJ(2,1));
t3076 = pkin(2) * t3104 - t3111 * t3110;
t3103 = sin(qJ(3,1));
t3168 = t3095 * t3103;
t3047 = pkin(3) * t3168 + t3076 * t3093;
t3109 = cos(qJ(3,1));
t3174 = t3093 * t3104;
t3197 = pkin(3) * t3109 ^ 2;
t3032 = 0.1e1 / (pkin(2) * t3168 + t3047 * t3109 + t3174 * t3197);
t3098 = legFrame(1,2);
t3085 = sin(t3098);
t3088 = cos(t3098);
t3073 = t3088 * g(1) - t3085 * g(2);
t3206 = t3032 * ((t3066 * t3085 - t3067 * t3088 + t3080) * t3110 + t3104 * (t3073 * t3092 + t3082));
t3202 = pkin(2) * t3099;
t3201 = pkin(2) * t3101;
t3200 = pkin(2) * t3103;
t3196 = pkin(3) * t3105;
t3195 = pkin(3) * t3107;
t3194 = pkin(3) * t3109;
t3165 = t3095 * t3106;
t3051 = t3092 * t3100 - t3094 * t3165;
t3077 = pkin(2) * t3106 + t3100 * t3111;
t3193 = (t3051 * t3196 + t3092 * t3074 - t3077 * t3173) * t3030;
t3163 = t3095 * t3108;
t3052 = t3092 * t3102 - t3094 * t3163;
t3078 = pkin(2) * t3108 + t3102 * t3111;
t3192 = (t3052 * t3195 + t3092 * t3075 - t3078 * t3173) * t3031;
t3161 = t3095 * t3110;
t3053 = t3092 * t3104 - t3094 * t3161;
t3079 = pkin(2) * t3110 + t3104 * t3111;
t3191 = (t3053 * t3194 + t3092 * t3076 - t3079 * t3173) * t3032;
t3171 = t3095 * t3100;
t3054 = t3092 * t3106 + t3094 * t3171;
t3148 = t3105 * t3093;
t3039 = -t3099 * t3054 - t3094 * t3148;
t3190 = t3030 * t3039;
t3068 = t3083 * g(1) + t3086 * g(2);
t3189 = t3030 * t3068;
t3188 = t3030 * t3083;
t3187 = t3030 * t3086;
t3169 = t3095 * t3102;
t3055 = t3092 * t3108 + t3094 * t3169;
t3146 = t3107 * t3093;
t3040 = -t3101 * t3055 - t3094 * t3146;
t3186 = t3031 * t3040;
t3069 = t3084 * g(1) + t3087 * g(2);
t3185 = t3031 * t3069;
t3184 = t3031 * t3084;
t3183 = t3031 * t3087;
t3167 = t3095 * t3104;
t3056 = t3092 * t3110 + t3094 * t3167;
t3144 = t3109 * t3093;
t3041 = -t3056 * t3103 - t3094 * t3144;
t3182 = t3032 * t3041;
t3070 = t3085 * g(1) + t3088 * g(2);
t3181 = t3032 * t3070;
t3180 = t3032 * t3085;
t3179 = t3032 * t3088;
t3177 = t3093 * t3094;
t3166 = t3095 * t3105;
t3164 = t3095 * t3107;
t3162 = t3095 * t3109;
t3160 = t3099 * t3093;
t3159 = t3099 * t3100;
t3158 = t3099 * t3106;
t3156 = t3101 * t3093;
t3155 = t3101 * t3102;
t3154 = t3101 * t3108;
t3152 = t3103 * t3093;
t3151 = t3103 * t3104;
t3150 = t3103 * t3110;
t3147 = t3105 * t3106;
t3145 = t3107 * t3108;
t3143 = t3109 * t3110;
t3142 = t3099 * t3208;
t3141 = t3101 * t3207;
t3140 = t3103 * t3206;
t3139 = t3105 * t3208;
t3138 = t3107 * t3207;
t3137 = t3109 * t3206;
t3057 = t3092 * t3165 + t3094 * t3100;
t3027 = t3057 * t3196 + t3074 * t3094 + t3077 * t3178;
t3136 = t3027 * t3188;
t3135 = t3027 * t3187;
t3058 = t3092 * t3163 + t3094 * t3102;
t3028 = t3058 * t3195 + t3075 * t3094 + t3078 * t3178;
t3134 = t3028 * t3184;
t3133 = t3028 * t3183;
t3059 = t3092 * t3161 + t3094 * t3104;
t3029 = t3059 * t3194 + t3076 * t3094 + t3079 * t3178;
t3132 = t3029 * t3180;
t3131 = t3029 * t3179;
t3048 = t3092 * t3171 - t3094 * t3106;
t3036 = t3099 * t3048 + t3092 * t3148;
t3130 = t3036 * t3188;
t3129 = t3036 * t3187;
t3049 = t3092 * t3169 - t3094 * t3108;
t3037 = t3101 * t3049 + t3092 * t3146;
t3128 = t3037 * t3184;
t3127 = t3037 * t3183;
t3050 = t3092 * t3167 - t3094 * t3110;
t3038 = t3103 * t3050 + t3092 * t3144;
t3126 = t3038 * t3180;
t3125 = t3038 * t3179;
t3124 = t3036 * t3142;
t3123 = t3037 * t3141;
t3122 = t3038 * t3140;
t3121 = t3036 * t3139;
t3120 = t3037 * t3138;
t3119 = t3038 * t3137;
t3118 = pkin(3) * t3160 - t3074 * t3095;
t3117 = pkin(3) * t3156 - t3075 * t3095;
t3116 = pkin(3) * t3152 - t3076 * t3095;
t3112 = 0.1e1 / pkin(3);
t3065 = t3104 * t3162 - t3152;
t3064 = t3102 * t3164 - t3156;
t3063 = t3100 * t3166 - t3160;
t3062 = t3095 * t3151 + t3144;
t3061 = t3095 * t3155 + t3146;
t3060 = t3095 * t3159 + t3148;
t3035 = -t3092 * t3079 + t3116 * t3094;
t3034 = -t3092 * t3078 + t3117 * t3094;
t3033 = -t3092 * t3077 + t3118 * t3094;
t3023 = -g(3) * t3050 + t3073 * t3056 + t3070 * t3174;
t3022 = -g(3) * t3049 + t3072 * t3055 + t3069 * t3175;
t3021 = -g(3) * t3048 + t3071 * t3054 + t3068 * t3176;
t3020 = -t3070 * t3093 * t3110 + g(3) * t3059 + t3073 * t3053;
t3019 = -t3069 * t3093 * t3108 + g(3) * t3058 + t3072 * t3052;
t3018 = -t3068 * t3093 * t3106 + g(3) * t3057 + t3071 * t3051;
t3011 = g(3) * (-t3060 * t3092 + t3094 * t3158) + t3071 * (t3060 * t3094 + t3092 * t3158) - t3068 * (-t3093 * t3159 + t3166);
t3010 = (-t3065 * t3092 + t3094 * t3143) * g(3) + (t3065 * t3094 + t3092 * t3143) * t3073 + t3070 * (t3104 * t3144 + t3168);
t3009 = (-t3064 * t3092 + t3094 * t3145) * g(3) + (t3064 * t3094 + t3092 * t3145) * t3072 + t3069 * (t3102 * t3146 + t3170);
t3008 = (-t3063 * t3092 + t3094 * t3147) * g(3) + (t3063 * t3094 + t3092 * t3147) * t3071 + t3068 * (t3100 * t3148 + t3172);
t3007 = (-t3062 * t3092 + t3094 * t3150) * g(3) + t3073 * (t3062 * t3094 + t3092 * t3150) - t3070 * (-t3093 * t3151 + t3162);
t3006 = (-t3061 * t3092 + t3094 * t3154) * g(3) + t3072 * (t3061 * t3094 + t3092 * t3154) - t3069 * (-t3093 * t3155 + t3164);
t1 = [(-((t3056 * t3088 + t3085 * t3174) * t3197 + (-t3035 * t3088 + t3047 * t3085) * t3109 + (t3095 * t3085 - t3088 * t3177) * t3200) * t3181 - ((t3055 * t3087 + t3084 * t3175) * t3198 + (-t3034 * t3087 + t3046 * t3084) * t3107 + (t3095 * t3084 - t3087 * t3177) * t3201) * t3185 - ((t3054 * t3086 + t3083 * t3176) * t3199 + (-t3033 * t3086 + t3045 * t3083) * t3105 + (t3095 * t3083 - t3086 * t3177) * t3202) * t3189) * MDP(1) + (-t3018 * t3129 - t3019 * t3127 - t3020 * t3125) * MDP(3) + (-t3021 * t3129 - t3022 * t3127 - t3023 * t3125) * MDP(4) + (-t3086 * t3121 - t3087 * t3120 - t3088 * t3119) * MDP(10) + (t3086 * t3124 + t3087 * t3123 + t3088 * t3122) * MDP(11) - g(1) * MDP(12) + ((-t3006 * t3133 - t3007 * t3131 - t3011 * t3135) * MDP(10) + (-t3008 * t3135 - t3009 * t3133 - t3010 * t3131) * MDP(11)) * t3112; (-(-(t3056 * t3085 - t3088 * t3174) * t3197 + (t3035 * t3085 + t3088 * t3047) * t3109 + (t3085 * t3177 + t3088 * t3095) * t3200) * t3181 - (-(t3055 * t3084 - t3087 * t3175) * t3198 + (t3034 * t3084 + t3087 * t3046) * t3107 + (t3084 * t3177 + t3087 * t3095) * t3201) * t3185 - (-(t3054 * t3083 - t3086 * t3176) * t3199 + (t3033 * t3083 + t3086 * t3045) * t3105 + (t3083 * t3177 + t3086 * t3095) * t3202) * t3189) * MDP(1) + (t3018 * t3130 + t3019 * t3128 + t3020 * t3126) * MDP(3) + (t3021 * t3130 + t3022 * t3128 + t3023 * t3126) * MDP(4) + (t3083 * t3121 + t3084 * t3120 + t3085 * t3119) * MDP(10) + (-t3083 * t3124 - t3084 * t3123 - t3085 * t3122) * MDP(11) - g(2) * MDP(12) + ((t3006 * t3134 + t3007 * t3132 + t3011 * t3136) * MDP(10) + (t3008 * t3136 + t3009 * t3134 + t3010 * t3132) * MDP(11)) * t3112; (-(-t3050 * t3197 + t3094 * t3079 * t3109 + (pkin(2) * t3152 + t3116 * t3109) * t3092) * t3181 - (-t3049 * t3198 + t3094 * t3078 * t3107 + (pkin(2) * t3156 + t3117 * t3107) * t3092) * t3185 - (-t3048 * t3199 + t3094 * t3077 * t3105 + (pkin(2) * t3160 + t3118 * t3105) * t3092) * t3189) * MDP(1) + (t3018 * t3190 + t3019 * t3186 + t3020 * t3182) * MDP(3) + (t3021 * t3190 + t3022 * t3186 + t3023 * t3182) * MDP(4) + (t3039 * t3139 + t3040 * t3138 + t3041 * t3137 + (t3006 * t3192 + t3007 * t3191 + t3011 * t3193) * t3112) * MDP(10) + (-t3039 * t3142 - t3040 * t3141 - t3041 * t3140 + (t3008 * t3193 + t3009 * t3192 + t3010 * t3191) * t3112) * MDP(11) - g(3) * MDP(12);];
taugX  = t1;
