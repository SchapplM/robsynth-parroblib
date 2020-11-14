% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:02:46
% EndTime: 2020-08-06 19:02:48
% DurationCPUTime: 1.70s
% Computational Cost: add. (927->140), mult. (1659->246), div. (132->6), fcn. (1641->18), ass. (0->120)
t3161 = MDP(10) - MDP(13);
t3097 = sin(qJ(1,3));
t3103 = cos(qJ(1,3));
t3074 = t3097 * g(1) - t3103 * g(2);
t3075 = t3103 * g(1) + t3097 * g(2);
t3093 = legFrame(3,3);
t3084 = sin(t3093);
t3087 = cos(t3093);
t3031 = t3074 * t3084 - t3075 * t3087;
t3096 = sin(qJ(2,3));
t3102 = cos(qJ(2,3));
t3009 = g(3) * t3102 + t3031 * t3096;
t3109 = 0.1e1 / qJ(3,3);
t3160 = t3009 * t3109;
t3099 = sin(qJ(1,2));
t3105 = cos(qJ(1,2));
t3076 = t3099 * g(1) - t3105 * g(2);
t3077 = t3105 * g(1) + t3099 * g(2);
t3094 = legFrame(2,3);
t3085 = sin(t3094);
t3088 = cos(t3094);
t3034 = t3076 * t3085 - t3077 * t3088;
t3098 = sin(qJ(2,2));
t3104 = cos(qJ(2,2));
t3013 = g(3) * t3104 + t3034 * t3098;
t3110 = 0.1e1 / qJ(3,2);
t3159 = t3013 * t3110;
t3095 = legFrame(1,3);
t3086 = sin(t3095);
t3089 = cos(t3095);
t3101 = sin(qJ(1,1));
t3107 = cos(qJ(1,1));
t3036 = (g(1) * t3101 - g(2) * t3107) * t3086 - (g(1) * t3107 + g(2) * t3101) * t3089;
t3100 = sin(qJ(2,1));
t3106 = cos(qJ(2,1));
t3017 = g(3) * t3106 + t3036 * t3100;
t3111 = 0.1e1 / qJ(3,1);
t3158 = t3017 * t3111;
t3156 = -MDP(12) + MDP(3);
t3010 = -g(3) * t3096 + t3031 * t3102;
t3014 = -g(3) * t3098 + t3034 * t3104;
t3018 = -g(3) * t3100 + t3036 * t3106;
t3052 = t3086 * g(1) - t3089 * g(2);
t3055 = t3089 * g(1) + t3086 * g(2);
t3021 = t3101 * t3052 - t3055 * t3107;
t3083 = t3100 * qJ(3,1);
t3073 = t3106 * pkin(1) + t3083;
t3120 = t3052 * t3107 + t3055 * t3101;
t3155 = t3156 * t3021 + (-MDP(14) * t3073 + t3161 * t3100 - MDP(2)) * t3120;
t3050 = t3084 * g(1) - t3087 * g(2);
t3053 = t3087 * g(1) + t3084 * g(2);
t3019 = t3097 * t3050 - t3053 * t3103;
t3029 = t3074 * t3087 + t3075 * t3084;
t3081 = t3096 * qJ(3,3);
t3071 = t3102 * pkin(1) + t3081;
t3153 = (MDP(10) * t3096 - MDP(2)) * (t3050 * t3103 + t3053 * t3097) - (MDP(13) * t3096 + MDP(14) * t3071) * t3029 + t3156 * t3019;
t3051 = t3085 * g(1) - t3088 * g(2);
t3054 = t3088 * g(1) + t3085 * g(2);
t3020 = t3099 * t3051 - t3054 * t3105;
t3032 = t3076 * t3088 + t3077 * t3085;
t3082 = t3098 * qJ(3,2);
t3072 = t3104 * pkin(1) + t3082;
t3152 = (MDP(10) * t3098 - MDP(2)) * (t3051 * t3105 + t3054 * t3099) - (MDP(13) * t3098 + MDP(14) * t3072) * t3032 + t3156 * t3020;
t3145 = t3097 * pkin(4);
t3144 = t3099 * pkin(4);
t3143 = t3101 * pkin(4);
t3142 = t3102 * qJ(3,3);
t3141 = t3104 * qJ(3,2);
t3140 = t3106 * qJ(3,1);
t3139 = MDP(9) + MDP(11);
t3108 = pkin(1) + pkin(2);
t3126 = t3108 * t3102;
t3065 = t3081 + t3126;
t3062 = 0.1e1 / t3065;
t3132 = t3062 * t3102;
t3125 = t3108 * t3104;
t3066 = t3082 + t3125;
t3063 = 0.1e1 / t3066;
t3131 = t3063 * t3104;
t3124 = t3108 * t3106;
t3067 = t3083 + t3124;
t3064 = 0.1e1 / t3067;
t3130 = t3064 * t3106;
t3129 = t3096 * t3109;
t3128 = t3098 * t3110;
t3127 = t3100 * t3111;
t2998 = -g(3) * t3071 - t3019 * (t3096 * pkin(1) - t3142);
t3114 = (MDP(14) * t2998 - t3161 * t3010) * t3109 * t3102;
t2999 = -g(3) * t3072 - t3020 * (t3098 * pkin(1) - t3141);
t3113 = (MDP(14) * t2999 - t3161 * t3014) * t3110 * t3104;
t3000 = -g(3) * t3073 - t3021 * (t3100 * pkin(1) - t3140);
t3112 = (MDP(14) * t3000 - t3161 * t3018) * t3111 * t3106;
t3092 = t3107 * pkin(4);
t3091 = t3105 * pkin(4);
t3090 = t3103 * pkin(4);
t3061 = t3107 * t3083 - t3143;
t3060 = t3101 * t3083 + t3092;
t3059 = t3105 * t3082 - t3144;
t3058 = t3099 * t3082 + t3091;
t3057 = t3103 * t3081 - t3145;
t3056 = t3097 * t3081 + t3090;
t3048 = t3086 * t3107 + t3089 * t3101;
t3047 = -t3086 * t3101 + t3089 * t3107;
t3046 = t3085 * t3105 + t3088 * t3099;
t3045 = -t3085 * t3099 + t3088 * t3105;
t3044 = t3084 * t3103 + t3087 * t3097;
t3043 = -t3084 * t3097 + t3087 * t3103;
t3042 = t3067 * t3107 - t3143;
t3041 = t3066 * t3105 - t3144;
t3040 = t3065 * t3103 - t3145;
t3039 = t3067 * t3101 + t3092;
t3038 = t3066 * t3099 + t3091;
t3037 = t3065 * t3097 + t3090;
t3006 = t3048 * t3124 + t3060 * t3089 + t3086 * t3061;
t3005 = t3047 * t3124 - t3086 * t3060 + t3061 * t3089;
t3004 = t3046 * t3125 + t3058 * t3088 + t3085 * t3059;
t3003 = t3045 * t3125 - t3085 * t3058 + t3059 * t3088;
t3002 = t3044 * t3126 + t3056 * t3087 + t3084 * t3057;
t3001 = t3043 * t3126 - t3084 * t3056 + t3057 * t3087;
t1 = [((-t3086 * t3039 + t3089 * t3042) * t3158 + (-t3085 * t3038 + t3088 * t3041) * t3159 + (-t3084 * t3037 + t3087 * t3040) * t3160) * MDP(14) - g(1) * MDP(15) + t3139 * ((-t3005 * t3158 - t3048 * t3120) * t3130 + (-t3003 * t3159 - t3032 * t3046) * t3131 + (-t3001 * t3160 - t3029 * t3044) * t3132) + (t3005 * t3112 + t3155 * t3048) * t3064 + (t3003 * t3113 + t3152 * t3046) * t3063 + (t3001 * t3114 + t3153 * t3044) * t3062; ((t3039 * t3089 + t3086 * t3042) * t3158 + (t3038 * t3088 + t3085 * t3041) * t3159 + (t3037 * t3087 + t3084 * t3040) * t3160) * MDP(14) - g(2) * MDP(15) + t3139 * ((-t3006 * t3158 + t3047 * t3120) * t3130 + (-t3004 * t3159 + t3032 * t3045) * t3131 + (-t3002 * t3160 + t3029 * t3043) * t3132) + (t3006 * t3112 - t3155 * t3047) * t3064 + (t3004 * t3113 - t3152 * t3045) * t3063 + (t3002 * t3114 - t3153 * t3043) * t3062; ((t3100 * t3000 + (t3108 * t3100 - t3140) * t3017) * t3111 + (t3098 * t2999 + (t3108 * t3098 - t3141) * t3013) * t3110 + (t3096 * t2998 + (t3108 * t3096 - t3142) * t3009) * t3109) * MDP(14) - g(3) * MDP(15) + t3139 * (-t3009 * t3129 - t3013 * t3128 - t3017 * t3127) - t3161 * (t3010 * t3129 + t3014 * t3128 + t3018 * t3127);];
taugX  = t1;
