% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:04:10
% EndTime: 2020-08-06 21:04:11
% DurationCPUTime: 1.42s
% Computational Cost: add. (702->161), mult. (975->267), div. (93->9), fcn. (909->35), ass. (0->116)
t3070 = MDP(3) - MDP(11);
t3014 = legFrame(3,3);
t2995 = sin(t3014);
t2998 = cos(t3014);
t2973 = g(1) * t2995 - g(2) * t2998;
t2976 = g(1) * t2998 + g(2) * t2995;
t3021 = sin(qJ(1,3));
t3027 = cos(qJ(1,3));
t2943 = t2973 * t3027 + t2976 * t3021;
t3017 = pkin(5) + qJ(3,3);
t3007 = -pkin(6) - t3017;
t3001 = 0.1e1 / t3007;
t3057 = t2943 * t3001;
t3015 = legFrame(2,3);
t2996 = sin(t3015);
t2999 = cos(t3015);
t2974 = g(1) * t2996 - g(2) * t2999;
t2977 = g(1) * t2999 + g(2) * t2996;
t3023 = sin(qJ(1,2));
t3029 = cos(qJ(1,2));
t2944 = t2974 * t3029 + t2977 * t3023;
t3018 = pkin(5) + qJ(3,2);
t3008 = -pkin(6) - t3018;
t3002 = 0.1e1 / t3008;
t3056 = t2944 * t3002;
t3016 = legFrame(1,3);
t2997 = sin(t3016);
t3000 = cos(t3016);
t2975 = g(1) * t2997 - g(2) * t3000;
t2978 = g(1) * t3000 + g(2) * t2997;
t3025 = sin(qJ(1,1));
t3031 = cos(qJ(1,1));
t2945 = t2975 * t3031 + t2978 * t3025;
t3019 = pkin(5) + qJ(3,1);
t3009 = -pkin(6) - t3019;
t3003 = 0.1e1 / t3009;
t3055 = t2945 * t3003;
t3010 = qJ(2,3) + pkin(7);
t3066 = pkin(3) * cos(t3010);
t3011 = qJ(2,2) + pkin(7);
t3065 = pkin(3) * cos(t3011);
t3012 = qJ(2,1) + pkin(7);
t3064 = pkin(3) * cos(t3012);
t3063 = sin(pkin(7)) * pkin(3);
t3062 = 0.2e1 * pkin(2) * pkin(3);
t3061 = 2 * pkin(1);
t2952 = (g(1) * t3027 + g(2) * t3021) * t2998 - (g(1) * t3021 - g(2) * t3027) * t2995;
t3026 = cos(qJ(2,3));
t3004 = t3026 * pkin(2);
t2979 = 0.1e1 / (t3004 + t3066);
t3020 = sin(qJ(2,3));
t3060 = (-g(3) * t3026 + t2952 * t3020) * t2979;
t2953 = (g(1) * t3029 + g(2) * t3023) * t2999 - (g(1) * t3023 - g(2) * t3029) * t2996;
t3028 = cos(qJ(2,2));
t3005 = t3028 * pkin(2);
t2980 = 0.1e1 / (t3005 + t3065);
t3022 = sin(qJ(2,2));
t3059 = (-g(3) * t3028 + t2953 * t3022) * t2980;
t2954 = (g(1) * t3031 + g(2) * t3025) * t3000 - (g(1) * t3025 - g(2) * t3031) * t2997;
t3030 = cos(qJ(2,1));
t3006 = t3030 * pkin(2);
t2981 = 0.1e1 / (t3006 + t3064);
t3024 = sin(qJ(2,1));
t3058 = (-g(3) * t3030 + t2954 * t3024) * t2981;
t2967 = -t2995 * t3021 + t2998 * t3027;
t3054 = t2967 * t3001;
t2968 = -t2996 * t3023 + t2999 * t3029;
t3053 = t2968 * t3002;
t2969 = -t2997 * t3025 + t3000 * t3031;
t3052 = t2969 * t3003;
t2970 = t2995 * t3027 + t2998 * t3021;
t3051 = t2970 * t3001;
t2971 = t2996 * t3029 + t2999 * t3023;
t3050 = t2971 * t3002;
t2972 = t2997 * t3031 + t3000 * t3025;
t3049 = t2972 * t3003;
t3048 = t3020 * t3057;
t3047 = t3026 * t3057;
t3046 = t3022 * t3056;
t3045 = t3028 * t3056;
t3044 = t3024 * t3055;
t3043 = t3030 * t3055;
t2988 = cos(pkin(7)) * pkin(3) + pkin(2);
t3042 = 0.1e1 / (t2988 * t3026 - t3020 * t3063) * (t2988 * t3020 + t3026 * t3063) * t3001;
t3041 = 0.1e1 / (t2988 * t3028 - t3022 * t3063) * (t2988 * t3022 + t3028 * t3063) * t3002;
t3040 = 0.1e1 / (t2988 * t3030 - t3024 * t3063) * (t2988 * t3024 + t3030 * t3063) * t3003;
t3039 = t2943 * t3042;
t3038 = t2944 * t3041;
t3037 = t2945 * t3040;
t2940 = t2973 * t3021 - t2976 * t3027;
t2941 = t2974 * t3023 - t2977 * t3029;
t2942 = t2975 * t3025 - t2978 * t3031;
t3036 = pkin(2) ^ 2;
t3035 = pkin(3) ^ 2;
t3034 = 0.2e1 * qJ(2,1);
t3033 = 0.2e1 * qJ(2,2);
t3032 = 0.2e1 * qJ(2,3);
t2991 = t3006 + pkin(1);
t2990 = t3005 + pkin(1);
t2989 = t3004 + pkin(1);
t2987 = t2990 * t3029;
t2986 = t2989 * t3027;
t2985 = t3031 * t2991;
t2984 = t3025 * t2991;
t2983 = t3023 * t2990;
t2982 = t3021 * t2989;
t2966 = -t3008 * t3023 + t2987;
t2965 = -t3007 * t3021 + t2986;
t2964 = -t3009 * t3025 + t2985;
t2963 = t3009 * t3031 + t2984;
t2962 = t3008 * t3029 + t2983;
t2961 = t3007 * t3027 + t2982;
t2936 = t2978 * (-t3019 * t3031 + t2984) + t2975 * (t3019 * t3025 + t2985);
t2935 = t2977 * (-t3018 * t3029 + t2983) + t2974 * (t3018 * t3023 + t2987);
t2934 = t2976 * (-t3017 * t3027 + t2982) + t2973 * (t3017 * t3021 + t2986);
t1 = [(-t2943 * t3054 - t2944 * t3053 - t2945 * t3052) * MDP(2) + (-t2967 * t3047 - t2968 * t3045 - t2969 * t3043) * MDP(9) + (t2967 * t3048 + t2968 * t3046 + t2969 * t3044) * MDP(10) + (-(t2969 * t2936 - (-t2963 * t2997 + t2964 * t3000 + t2969 * t3064) * t2945) * t3003 - (t2968 * t2935 - (-t2962 * t2996 + t2966 * t2999 + t2968 * t3065) * t2944) * t3002 - (t2967 * t2934 - (-t2961 * t2995 + t2965 * t2998 + t2967 * t3066) * t2943) * t3001) * MDP(12) - g(1) * MDP(13) + t3070 * (t2940 * t3054 + t2941 * t3053 + t2942 * t3052); (-t2943 * t3051 - t2944 * t3050 - t2945 * t3049) * MDP(2) + (-t2970 * t3047 - t2971 * t3045 - t2972 * t3043) * MDP(9) + (t2970 * t3048 + t2971 * t3046 + t2972 * t3044) * MDP(10) + (-(t2972 * t2936 - (t2963 * t3000 + t2964 * t2997 + t2972 * t3064) * t2945) * t3003 - (t2971 * t2935 - (t2962 * t2999 + t2966 * t2996 + t2971 * t3065) * t2944) * t3002 - (t2970 * t2934 - (t2961 * t2998 + t2965 * t2995 + t2970 * t3066) * t2943) * t3001) * MDP(12) - g(2) * MDP(13) + t3070 * (t2940 * t3051 + t2941 * t3050 + t2942 * t3049); (-t3037 - t3038 - t3039) * MDP(2) + (-t3026 * t3039 - t3028 * t3038 - t3030 * t3037 + t3058 + t3059 + t3060) * MDP(9) + (t3024 * t3037 + t2981 * (g(3) * t3024 + t2954 * t3030) + t3022 * t3038 + t2980 * (g(3) * t3022 + t2953 * t3028) + t3020 * t3039 + t2979 * (g(3) * t3020 + t2952 * t3026)) * MDP(10) + (-t2936 * t3040 + pkin(2) * t3058 + (sin(t3034 + pkin(7)) * t3062 + t3036 * sin(t3034) + t3035 * sin(0.2e1 * t3012) + (sin(t3012) * pkin(3) + pkin(2) * t3024) * t3061) * t2981 * t3055 / 0.2e1 - t2935 * t3041 + pkin(2) * t3059 + (sin(t3033 + pkin(7)) * t3062 + t3036 * sin(t3033) + t3035 * sin(0.2e1 * t3011) + (sin(t3011) * pkin(3) + pkin(2) * t3022) * t3061) * t2980 * t3056 / 0.2e1 - t2934 * t3042 + pkin(2) * t3060 + (sin(t3032 + pkin(7)) * t3062 + t3036 * sin(t3032) + t3035 * sin(0.2e1 * t3010) + (sin(t3010) * pkin(3) + pkin(2) * t3020) * t3061) * t2979 * t3057 / 0.2e1) * MDP(12) - g(3) * MDP(13) + t3070 * (t2940 * t3042 + t2941 * t3041 + t2942 * t3040);];
taugX  = t1;
