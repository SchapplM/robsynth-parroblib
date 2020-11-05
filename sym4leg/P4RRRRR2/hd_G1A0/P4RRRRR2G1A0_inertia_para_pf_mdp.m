% Calculate minimal parameter regressor of inertia matrix for parallel robot
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [17x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4RRRRR2G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [4x4]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P4RRRRR2G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(17,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_mdp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:26:06
% EndTime: 2020-08-07 17:26:35
% DurationCPUTime: 30.34s
% Computational Cost: add. (21066->1194), mult. (12609->1867), div. (8008->24), fcn. (14292->66), ass. (0->649)
t3234 = sin(qJ(3,1));
t3240 = cos(qJ(3,1));
t3712 = MDP(12) * t3240 - MDP(13) * t3234;
t3232 = sin(qJ(3,2));
t3238 = cos(qJ(3,2));
t3711 = MDP(12) * t3238 - MDP(13) * t3232;
t3230 = sin(qJ(3,3));
t3236 = cos(qJ(3,3));
t3710 = MDP(12) * t3236 - MDP(13) * t3230;
t3226 = sin(qJ(3,4));
t3228 = cos(qJ(3,4));
t3709 = MDP(12) * t3228 - MDP(13) * t3226;
t3708 = 2 * pkin(1);
t3707 = 2 * MDP(8);
t3202 = 0.1e1 / t3228;
t3217 = 0.1e1 / t3236;
t3220 = 0.1e1 / t3238;
t3223 = 0.1e1 / t3240;
t3203 = 0.1e1 / t3228 ^ 2;
t3218 = 0.1e1 / t3236 ^ 2;
t3221 = 0.1e1 / t3238 ^ 2;
t3224 = 0.1e1 / t3240 ^ 2;
t3253 = 1 / pkin(1);
t3227 = sin(qJ(2,4));
t3200 = 0.1e1 / t3227;
t3597 = t3202 * t3226;
t3482 = t3200 * t3597;
t3147 = t3253 * t3482;
t3229 = cos(qJ(2,4));
t3148 = pkin(1) * t3229 + pkin(2) * t3228;
t3621 = t3148 * t3203;
t3507 = t3226 * t3621;
t3400 = t3200 * t3507;
t3251 = 0.1e1 / pkin(2);
t3565 = t3251 * t3253;
t3274 = t3400 * t3565;
t3094 = t3147 - t3274 / 0.2e1;
t3576 = t3226 * t3229;
t3698 = t3251 / 0.2e1;
t3706 = -0.2e1 * t3094 * t3576 - 0.2e1 * t3227 * t3698;
t3231 = sin(qJ(2,3));
t3207 = 0.1e1 / t3231;
t3583 = t3217 * t3230;
t3472 = t3207 * t3583;
t3149 = t3253 * t3472;
t3237 = cos(qJ(2,3));
t3154 = pkin(1) * t3237 + pkin(2) * t3236;
t3619 = t3154 * t3218;
t3505 = t3230 * t3619;
t3394 = t3207 * t3505;
t3273 = t3394 * t3565;
t3099 = t3149 - t3273 / 0.2e1;
t3573 = t3230 * t3237;
t3705 = -0.2e1 * t3099 * t3573 - 0.2e1 * t3231 * t3698;
t3233 = sin(qJ(2,2));
t3211 = 0.1e1 / t3233;
t3581 = t3220 * t3232;
t3466 = t3211 * t3581;
t3150 = t3253 * t3466;
t3239 = cos(qJ(2,2));
t3155 = pkin(1) * t3239 + pkin(2) * t3238;
t3617 = t3155 * t3221;
t3503 = t3232 * t3617;
t3388 = t3211 * t3503;
t3272 = t3388 * t3565;
t3100 = t3150 - t3272 / 0.2e1;
t3571 = t3232 * t3239;
t3704 = -0.2e1 * t3100 * t3571 - 0.2e1 * t3233 * t3698;
t3235 = sin(qJ(2,1));
t3215 = 0.1e1 / t3235;
t3579 = t3223 * t3234;
t3460 = t3215 * t3579;
t3151 = t3253 * t3460;
t3241 = cos(qJ(2,1));
t3156 = pkin(1) * t3241 + pkin(2) * t3240;
t3615 = t3156 * t3224;
t3501 = t3234 * t3615;
t3382 = t3215 * t3501;
t3271 = t3382 * t3565;
t3101 = t3151 - t3271 / 0.2e1;
t3569 = t3234 * t3241;
t3703 = -0.2e1 * t3101 * t3569 - 0.2e1 * t3235 * t3698;
t3441 = qJ(1,4) + legFrame(4,3);
t3188 = qJ(2,4) + t3441;
t3242 = xP(4);
t3157 = -t3242 + t3188;
t3243 = koppelP(4,2);
t3247 = koppelP(4,1);
t3114 = -t3243 * cos(t3157) + t3247 * sin(t3157);
t3702 = 0.2e1 * t3114;
t3442 = qJ(1,3) + legFrame(3,3);
t3192 = qJ(2,3) + t3442;
t3165 = -t3242 + t3192;
t3244 = koppelP(3,2);
t3248 = koppelP(3,1);
t3115 = -t3244 * cos(t3165) + t3248 * sin(t3165);
t3701 = 0.2e1 * t3115;
t3443 = qJ(1,2) + legFrame(2,3);
t3193 = qJ(2,2) + t3443;
t3166 = -t3242 + t3193;
t3245 = koppelP(2,2);
t3249 = koppelP(2,1);
t3116 = -t3245 * cos(t3166) + t3249 * sin(t3166);
t3700 = 0.2e1 * t3116;
t3444 = qJ(1,1) + legFrame(1,3);
t3194 = qJ(2,1) + t3444;
t3167 = -t3242 + t3194;
t3246 = koppelP(1,2);
t3250 = koppelP(1,1);
t3117 = -t3246 * cos(t3167) + t3250 * sin(t3167);
t3699 = 0.2e1 * t3117;
t3697 = pkin(1) * sin(t3441);
t3696 = pkin(1) * sin(t3442);
t3695 = pkin(1) * sin(t3443);
t3694 = pkin(1) * sin(t3444);
t3693 = MDP(1) * t3253;
t3692 = MDP(1) / pkin(1) ^ 2;
t3691 = MDP(4) * t3253;
t3690 = MDP(7) * t3253;
t3689 = MDP(9) * t3251;
t3688 = MDP(9) * t3253;
t3687 = MDP(10) * t3251;
t3686 = MDP(10) * t3253;
t3196 = sin(t3242);
t3197 = cos(t3242);
t3126 = t3196 * t3247 + t3197 * t3243;
t3172 = cos(t3188);
t3598 = t3200 * t3253;
t3135 = t3172 * t3598;
t3118 = t3126 * t3135;
t3130 = -t3196 * t3243 + t3197 * t3247;
t3171 = sin(t3188);
t3134 = t3171 * t3598;
t3119 = t3130 * t3134;
t3173 = qJ(3,4) + t3188;
t3174 = -qJ(3,4) + t3188;
t3098 = cos(t3441) * t3708 + (cos(t3173) + cos(t3174)) * pkin(2);
t3143 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t3518 = t3143 * t3565;
t3084 = t3098 * t3518;
t3044 = t3126 * t3084;
t3152 = sin(t3173);
t3153 = sin(t3174);
t3097 = -0.2e1 * t3697 + (-t3152 - t3153) * pkin(2);
t3083 = t3097 * t3518;
t3045 = t3130 * t3083;
t3563 = t3044 + t3045;
t3024 = -0.2e1 * t3118 + 0.2e1 * t3119 + t3563;
t3677 = t3024 * t3114;
t3131 = -t3196 * t3244 + t3197 * t3248;
t3175 = sin(t3192);
t3592 = t3207 * t3253;
t3137 = t3175 * t3592;
t3120 = t3131 * t3137;
t3127 = t3196 * t3248 + t3197 * t3244;
t3178 = cos(t3192);
t3140 = t3178 * t3592;
t3123 = t3127 * t3140;
t3181 = qJ(3,3) + t3192;
t3182 = -qJ(3,3) + t3192;
t3111 = cos(t3442) * t3708 + (cos(t3181) + cos(t3182)) * pkin(2);
t3144 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t3515 = t3144 * t3565;
t3091 = t3111 * t3515;
t3050 = t3127 * t3091;
t3158 = sin(t3181);
t3159 = sin(t3182);
t3108 = -0.2e1 * t3696 + (-t3158 - t3159) * pkin(2);
t3088 = t3108 * t3515;
t3053 = t3131 * t3088;
t3562 = t3050 + t3053;
t3029 = -0.2e1 * t3123 + 0.2e1 * t3120 + t3562;
t3676 = t3029 * t3115;
t3132 = -t3196 * t3245 + t3197 * t3249;
t3176 = sin(t3193);
t3588 = t3211 * t3253;
t3138 = t3176 * t3588;
t3121 = t3132 * t3138;
t3128 = t3196 * t3249 + t3197 * t3245;
t3179 = cos(t3193);
t3141 = t3179 * t3588;
t3124 = t3128 * t3141;
t3183 = qJ(3,2) + t3193;
t3184 = -qJ(3,2) + t3193;
t3112 = cos(t3443) * t3708 + (cos(t3183) + cos(t3184)) * pkin(2);
t3145 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t3512 = t3145 * t3565;
t3092 = t3112 * t3512;
t3051 = t3128 * t3092;
t3160 = sin(t3183);
t3161 = sin(t3184);
t3109 = -0.2e1 * t3695 + (-t3160 - t3161) * pkin(2);
t3089 = t3109 * t3512;
t3054 = t3132 * t3089;
t3561 = t3051 + t3054;
t3031 = -0.2e1 * t3124 + 0.2e1 * t3121 + t3561;
t3675 = t3031 * t3116;
t3133 = -t3196 * t3246 + t3197 * t3250;
t3177 = sin(t3194);
t3584 = t3215 * t3253;
t3139 = t3177 * t3584;
t3122 = t3133 * t3139;
t3129 = t3196 * t3250 + t3197 * t3246;
t3180 = cos(t3194);
t3142 = t3180 * t3584;
t3125 = t3129 * t3142;
t3185 = qJ(3,1) + t3194;
t3186 = -qJ(3,1) + t3194;
t3113 = cos(t3444) * t3708 + (cos(t3185) + cos(t3186)) * pkin(2);
t3146 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t3509 = t3146 * t3565;
t3093 = t3113 * t3509;
t3052 = t3129 * t3093;
t3162 = sin(t3185);
t3163 = sin(t3186);
t3110 = -0.2e1 * t3694 + (-t3162 - t3163) * pkin(2);
t3090 = t3110 * t3509;
t3055 = t3133 * t3090;
t3560 = t3052 + t3055;
t3033 = -0.2e1 * t3125 + 0.2e1 * t3122 + t3560;
t3674 = t3033 * t3117;
t3035 = (t3126 * t3098 - 0.2e1 * t3130 * (t3697 + (t3153 / 0.2e1 + t3152 / 0.2e1) * pkin(2))) * t3518;
t3673 = t3035 * t3171;
t3672 = t3035 * t3172;
t3199 = t3226 ^ 2;
t3671 = t3035 * t3199;
t3036 = (t3127 * t3111 - 0.2e1 * t3131 * (t3696 + (t3159 / 0.2e1 + t3158 / 0.2e1) * pkin(2))) * t3515;
t3670 = t3036 * t3175;
t3669 = t3036 * t3178;
t3206 = t3230 ^ 2;
t3668 = t3036 * t3206;
t3037 = (t3128 * t3112 - 0.2e1 * t3132 * (t3695 + (t3161 / 0.2e1 + t3160 / 0.2e1) * pkin(2))) * t3512;
t3667 = t3037 * t3176;
t3666 = t3037 * t3179;
t3210 = t3232 ^ 2;
t3665 = t3037 * t3210;
t3038 = (t3129 * t3113 - 0.2e1 * t3133 * (t3694 + (t3163 / 0.2e1 + t3162 / 0.2e1) * pkin(2))) * t3509;
t3664 = t3038 * t3177;
t3663 = t3038 * t3180;
t3214 = t3234 ^ 2;
t3662 = t3038 * t3214;
t3074 = -t3118 + t3119;
t3661 = t3074 * t3227;
t3660 = t3074 * t3229;
t3075 = -t3123 + t3120;
t3659 = t3075 * t3231;
t3658 = t3075 * t3237;
t3076 = -t3124 + t3121;
t3657 = t3076 * t3233;
t3656 = t3076 * t3239;
t3077 = -t3125 + t3122;
t3655 = t3077 * t3235;
t3654 = t3077 * t3241;
t3095 = 0.2e1 * t3147 - t3274;
t3653 = t3095 * t3114;
t3652 = t3097 * t3143;
t3651 = t3098 * t3143;
t3102 = 0.2e1 * t3149 - t3273;
t3650 = t3102 * t3115;
t3104 = 0.2e1 * t3150 - t3272;
t3649 = t3104 * t3116;
t3106 = 0.2e1 * t3151 - t3271;
t3648 = t3106 * t3117;
t3647 = t3108 * t3144;
t3646 = t3109 * t3145;
t3645 = t3110 * t3146;
t3644 = t3111 * t3144;
t3643 = t3112 * t3145;
t3642 = t3113 * t3146;
t3641 = t3114 * t3200;
t3201 = 0.1e1 / t3227 ^ 2;
t3640 = t3114 * t3201;
t3639 = t3115 * t3207;
t3208 = 0.1e1 / t3231 ^ 2;
t3638 = t3115 * t3208;
t3637 = t3116 * t3211;
t3212 = 0.1e1 / t3233 ^ 2;
t3636 = t3116 * t3212;
t3635 = t3117 * t3215;
t3216 = 0.1e1 / t3235 ^ 2;
t3634 = t3117 * t3216;
t3633 = t3143 * t3171;
t3632 = t3143 * t3172;
t3631 = t3143 * t3199;
t3630 = t3144 * t3175;
t3629 = t3144 * t3178;
t3628 = t3144 * t3206;
t3627 = t3145 * t3176;
t3626 = t3145 * t3179;
t3625 = t3145 * t3210;
t3624 = t3146 * t3177;
t3623 = t3146 * t3180;
t3622 = t3146 * t3214;
t3620 = t3148 * t3229;
t3618 = t3154 * t3237;
t3616 = t3155 * t3239;
t3614 = t3156 * t3241;
t3613 = t3171 * t3200;
t3612 = t3172 * t3200;
t3611 = t3172 * t3201;
t3610 = t3175 * t3207;
t3609 = t3176 * t3211;
t3608 = t3177 * t3215;
t3607 = t3178 * t3207;
t3606 = t3178 * t3208;
t3605 = t3179 * t3211;
t3604 = t3179 * t3212;
t3603 = t3180 * t3215;
t3602 = t3180 * t3216;
t3601 = t3199 * t3200;
t3600 = t3200 * t3202;
t3599 = t3200 * t3229;
t3596 = t3203 * t3199;
t3595 = t3206 * t3207;
t3594 = t3207 * t3217;
t3593 = t3207 * t3237;
t3591 = t3210 * t3211;
t3590 = t3211 * t3220;
t3589 = t3211 * t3239;
t3587 = t3214 * t3215;
t3586 = t3215 * t3223;
t3585 = t3215 * t3241;
t3582 = t3218 * t3206;
t3580 = t3221 * t3210;
t3578 = t3224 * t3214;
t3577 = t3226 * t3228;
t3575 = t3228 * t3229;
t3574 = t3230 * t3236;
t3572 = t3232 * t3238;
t3570 = t3234 * t3240;
t3568 = t3236 * t3237;
t3567 = t3238 * t3239;
t3566 = t3240 * t3241;
t3252 = 0.1e1 / pkin(2) ^ 2;
t3564 = t3252 * t3253;
t3559 = 0.2e1 * MDP(12);
t3558 = 0.2e1 * MDP(13);
t3557 = t3253 * t3707;
t3556 = t3200 * t3706;
t3555 = t3207 * t3705;
t3554 = t3211 * t3704;
t3553 = t3215 * t3703;
t3552 = t3035 * t3577;
t3551 = t3036 * t3574;
t3550 = t3037 * t3572;
t3549 = t3038 * t3570;
t3548 = t3143 * t3661;
t3547 = t3143 * t3660;
t3546 = t3074 * t3576;
t3545 = t3074 * t3575;
t3544 = t3144 * t3659;
t3543 = t3144 * t3658;
t3542 = t3075 * t3573;
t3541 = t3075 * t3568;
t3540 = t3145 * t3657;
t3539 = t3145 * t3656;
t3538 = t3076 * t3571;
t3537 = t3076 * t3567;
t3536 = t3146 * t3655;
t3535 = t3146 * t3654;
t3534 = t3077 * t3569;
t3533 = t3077 * t3566;
t3532 = t3097 * t3631;
t3531 = t3098 * t3631;
t3530 = t3108 * t3628;
t3529 = t3109 * t3625;
t3528 = t3110 * t3622;
t3527 = t3111 * t3628;
t3526 = t3112 * t3625;
t3525 = t3113 * t3622;
t3524 = t3114 * t3601;
t3523 = t3115 * t3595;
t3522 = t3116 * t3591;
t3521 = t3117 * t3587;
t3520 = t3143 * t3597;
t3519 = t3143 * t3577;
t3517 = t3144 * t3583;
t3516 = t3144 * t3574;
t3514 = t3145 * t3581;
t3513 = t3145 * t3572;
t3511 = t3146 * t3579;
t3510 = t3146 * t3570;
t3508 = t3201 * t3620;
t3506 = t3208 * t3618;
t3504 = t3212 * t3616;
t3502 = t3216 * t3614;
t3500 = t3171 * t3601;
t3499 = t3171 * t3599;
t3498 = t3172 * t3601;
t3497 = t3172 * t3599;
t3496 = t3175 * t3595;
t3495 = t3175 * t3593;
t3494 = t3176 * t3591;
t3493 = t3176 * t3589;
t3492 = t3177 * t3587;
t3491 = t3177 * t3585;
t3490 = t3178 * t3595;
t3489 = t3178 * t3593;
t3488 = t3179 * t3591;
t3487 = t3179 * t3589;
t3486 = t3180 * t3587;
t3485 = t3180 * t3585;
t3198 = t3226 * t3199;
t3484 = t3198 * t3600;
t3483 = t3199 * t3600;
t3481 = t3200 * t3577;
t3480 = t3200 * t3576;
t3479 = t3200 * t3575;
t3096 = t3147 - t3274;
t3478 = t3096 * t3597;
t3477 = t3171 * t3597;
t3476 = t3172 * t3597;
t3475 = t3227 * t3597;
t3205 = t3230 * t3206;
t3474 = t3205 * t3594;
t3473 = t3206 * t3594;
t3471 = t3207 * t3574;
t3470 = t3207 * t3573;
t3469 = t3207 * t3568;
t3209 = t3232 * t3210;
t3468 = t3209 * t3590;
t3467 = t3210 * t3590;
t3465 = t3211 * t3572;
t3464 = t3211 * t3571;
t3463 = t3211 * t3567;
t3213 = t3234 * t3214;
t3462 = t3213 * t3586;
t3461 = t3214 * t3586;
t3459 = t3215 * t3570;
t3458 = t3215 * t3569;
t3457 = t3215 * t3566;
t3103 = t3149 - t3273;
t3456 = t3103 * t3583;
t3455 = t3175 * t3583;
t3454 = t3178 * t3583;
t3453 = t3231 * t3583;
t3105 = t3150 - t3272;
t3452 = t3105 * t3581;
t3451 = t3176 * t3581;
t3450 = t3179 * t3581;
t3449 = t3233 * t3581;
t3107 = t3151 - t3271;
t3448 = t3107 * t3579;
t3447 = t3177 * t3579;
t3446 = t3180 * t3579;
t3445 = t3235 * t3579;
t3440 = t3197 * MDP(15) - t3196 * MDP(16);
t3023 = t3044 / 0.2e1 + t3045 / 0.2e1 + t3074;
t3439 = t3023 * t3479;
t3026 = t3050 / 0.2e1 + t3053 / 0.2e1 + t3075;
t3438 = t3026 * t3469;
t3027 = t3051 / 0.2e1 + t3054 / 0.2e1 + t3076;
t3437 = t3027 * t3463;
t3028 = t3052 / 0.2e1 + t3055 / 0.2e1 + t3077;
t3436 = t3028 * t3457;
t3435 = t3143 * t3546;
t3434 = t3143 * t3545;
t3433 = t3144 * t3542;
t3432 = t3144 * t3541;
t3431 = t3145 * t3538;
t3430 = t3145 * t3537;
t3429 = t3146 * t3534;
t3428 = t3146 * t3533;
t3427 = t3097 * t3519;
t3426 = t3098 * t3519;
t3425 = t3108 * t3516;
t3424 = t3109 * t3513;
t3423 = t3110 * t3510;
t3422 = t3111 * t3516;
t3421 = t3112 * t3513;
t3420 = t3113 * t3510;
t3419 = t3114 * t3481;
t3418 = t3115 * t3471;
t3417 = t3116 * t3465;
t3416 = t3117 * t3459;
t3415 = t3143 * t3499;
t3414 = t3143 * t3497;
t3413 = t3143 * t3480;
t3412 = t3144 * t3495;
t3411 = t3144 * t3489;
t3410 = t3144 * t3470;
t3409 = t3145 * t3493;
t3408 = t3145 * t3487;
t3407 = t3145 * t3464;
t3406 = t3146 * t3491;
t3405 = t3146 * t3485;
t3404 = t3146 * t3458;
t3403 = t3198 * t3200 * t3621;
t3402 = t3148 * t3483;
t3401 = t3596 * t3620;
t3204 = t3202 * t3203;
t3399 = t3204 * t3508;
t3398 = t3226 * t3508;
t3397 = t3205 * t3207 * t3619;
t3396 = t3154 * t3473;
t3395 = t3582 * t3618;
t3219 = t3217 * t3218;
t3393 = t3219 * t3506;
t3392 = t3230 * t3506;
t3391 = t3209 * t3211 * t3617;
t3390 = t3155 * t3467;
t3389 = t3580 * t3616;
t3222 = t3220 * t3221;
t3387 = t3222 * t3504;
t3386 = t3232 * t3504;
t3385 = t3213 * t3215 * t3615;
t3384 = t3156 * t3461;
t3383 = t3578 * t3614;
t3225 = t3223 * t3224;
t3381 = t3225 * t3502;
t3380 = t3234 * t3502;
t3379 = t3171 * t3481;
t3378 = t3171 * t3479;
t3377 = t3172 * t3481;
t3376 = t3172 * t3479;
t3375 = t3175 * t3471;
t3374 = t3175 * t3469;
t3373 = t3176 * t3465;
t3372 = t3176 * t3463;
t3371 = t3177 * t3459;
t3370 = t3177 * t3457;
t3369 = t3178 * t3471;
t3368 = t3178 * t3469;
t3367 = t3179 * t3465;
t3366 = t3179 * t3463;
t3365 = t3180 * t3459;
t3364 = t3180 * t3457;
t3363 = t3229 * t3483;
t3362 = t3202 * t3480;
t3361 = t3237 * t3473;
t3360 = t3217 * t3470;
t3359 = t3239 * t3467;
t3358 = t3220 * t3464;
t3357 = t3241 * t3461;
t3356 = t3223 * t3458;
t3355 = t3023 * t3480;
t3042 = t3134 + t3083 / 0.2e1;
t3354 = t3042 * t3480;
t3043 = t3135 - t3084 / 0.2e1;
t3353 = t3043 * t3480;
t3352 = t3026 * t3470;
t3056 = t3137 + t3088 / 0.2e1;
t3351 = t3056 * t3470;
t3059 = t3140 - t3091 / 0.2e1;
t3350 = t3059 * t3470;
t3349 = t3027 * t3464;
t3057 = t3138 + t3089 / 0.2e1;
t3348 = t3057 * t3464;
t3060 = t3141 - t3092 / 0.2e1;
t3347 = t3060 * t3464;
t3346 = t3028 * t3458;
t3058 = t3139 + t3090 / 0.2e1;
t3345 = t3058 * t3458;
t3061 = t3142 - t3093 / 0.2e1;
t3344 = t3061 * t3458;
t3343 = 0.2e1 * t3439;
t3342 = 0.2e1 * t3438;
t3341 = 0.2e1 * t3437;
t3340 = 0.2e1 * t3436;
t3339 = -0.2e1 * t3363;
t3338 = -0.2e1 * t3361;
t3337 = -0.2e1 * t3359;
t3336 = -0.2e1 * t3357;
t3335 = -0.2e1 * t3355;
t3334 = -0.2e1 * t3352;
t3333 = -0.2e1 * t3349;
t3332 = -0.2e1 * t3346;
t3331 = t3074 * t3148 * t3480;
t3330 = t3075 * t3154 * t3470;
t3329 = t3076 * t3155 * t3464;
t3328 = t3077 * t3156 * t3458;
t3327 = t3171 * t3413;
t3326 = t3143 * t3378;
t3325 = t3172 * t3413;
t3324 = t3143 * t3376;
t3323 = t3143 * t3363;
t3322 = t3143 * t3362;
t3321 = t3175 * t3410;
t3320 = t3144 * t3374;
t3319 = t3178 * t3410;
t3318 = t3144 * t3368;
t3317 = t3144 * t3361;
t3316 = t3144 * t3360;
t3315 = t3176 * t3407;
t3314 = t3145 * t3372;
t3313 = t3179 * t3407;
t3312 = t3145 * t3366;
t3311 = t3145 * t3359;
t3310 = t3145 * t3358;
t3309 = t3177 * t3404;
t3308 = t3146 * t3370;
t3307 = t3180 * t3404;
t3306 = t3146 * t3364;
t3305 = t3146 * t3357;
t3304 = t3146 * t3356;
t3303 = t3201 * t3401;
t3302 = t3202 * t3398;
t3301 = t3203 * t3398;
t3300 = t3208 * t3395;
t3299 = t3217 * t3392;
t3298 = t3218 * t3392;
t3297 = t3212 * t3389;
t3296 = t3220 * t3386;
t3295 = t3221 * t3386;
t3294 = t3216 * t3383;
t3293 = t3223 * t3380;
t3292 = t3224 * t3380;
t3291 = -t3196 * MDP(15) - t3197 * MDP(16);
t3046 = 0.2e1 * t3134 + t3083;
t3288 = t3046 * t3114 + t3673;
t3048 = 0.2e1 * t3135 - t3084;
t3287 = t3048 * t3114 + t3672;
t3062 = 0.2e1 * t3137 + t3088;
t3284 = t3062 * t3115 + t3670;
t3068 = 0.2e1 * t3140 - t3091;
t3283 = t3068 * t3115 + t3669;
t3064 = 0.2e1 * t3138 + t3089;
t3280 = t3064 * t3116 + t3667;
t3070 = 0.2e1 * t3141 - t3092;
t3279 = t3070 * t3116 + t3666;
t3066 = 0.2e1 * t3139 + t3090;
t3276 = t3066 * t3117 + t3664;
t3072 = 0.2e1 * t3142 - t3093;
t3275 = t3072 * t3117 + t3663;
t3270 = t3035 * t3597 + t3036 * t3583 + t3037 * t3581 + t3038 * t3579;
t3269 = t3097 * t3520 + t3108 * t3517 + t3109 * t3514 + t3110 * t3511;
t3268 = t3098 * t3520 + t3111 * t3517 + t3112 * t3514 + t3113 * t3511;
t3267 = t3148 * t3204 * t3601 + t3154 * t3219 * t3595 + t3155 * t3222 * t3591 + t3156 * t3225 * t3587;
t3136 = (t3196 ^ 2 + t3197 ^ 2) * MDP(17);
t3081 = -0.2e1 * t3101 * t3566 + t3251 * t3445;
t3080 = -0.2e1 * t3100 * t3567 + t3251 * t3449;
t3079 = -0.2e1 * t3099 * t3568 + t3251 * t3453;
t3078 = -0.2e1 * t3094 * t3575 + t3251 * t3475;
t3073 = t3142 - t3093;
t3071 = t3141 - t3092;
t3069 = t3140 - t3091;
t3067 = t3139 + t3090;
t3065 = t3138 + t3089;
t3063 = t3137 + t3088;
t3049 = t3135 - t3084;
t3047 = t3134 + t3083;
t3041 = (t3201 * t3476 + t3208 * t3454 + t3212 * t3450 + t3216 * t3446) * t3692;
t3040 = (t3201 * t3477 + t3208 * t3455 + t3212 * t3451 + t3216 * t3447) * t3692;
t3039 = (t3171 * t3611 + t3175 * t3606 + t3176 * t3604 + t3177 * t3602) * t3692;
t3034 = t3077 + t3560;
t3032 = t3076 + t3561;
t3030 = t3075 + t3562;
t3025 = t3074 + t3563;
t1 = [(t3048 * t3497 + t3068 * t3489 + t3070 * t3487 + t3072 * t3485) * MDP(5) + (-t3048 * t3172 - t3068 * t3178 - t3070 * t3179 - t3072 * t3180) * MDP(6) + (t3043 * t3376 + t3059 * t3368 + t3060 * t3366 + t3061 * t3364) * t3559 + (-t3172 * t3353 - t3178 * t3350 - t3179 * t3347 - t3180 * t3344) * t3558 + t3136 + (t3172 ^ 2 * t3201 + t3178 ^ 2 * t3208 + t3179 ^ 2 * t3212 + t3180 ^ 2 * t3216) * t3692 + ((t3049 * t3612 + t3069 * t3607 + t3071 * t3605 + t3073 * t3603) * MDP(4) + (t3049 * t3498 + t3069 * t3490 + t3071 * t3488 + t3073 * t3486) * MDP(7) + (t3049 * t3377 + t3069 * t3369 + t3071 * t3367 + t3073 * t3365) * t3707 + ((-t3049 * t3651 - t3069 * t3644 - t3071 * t3643 - t3073 * t3642) * MDP(4) + (-t3098 * t3414 - t3111 * t3411 - t3112 * t3408 - t3113 * t3405) * MDP(5) + (t3098 * t3632 + t3111 * t3629 + t3112 * t3626 + t3113 * t3623) * MDP(6) + (-t3049 * t3531 - t3069 * t3527 - t3071 * t3526 - t3073 * t3525) * MDP(7) + (-t3049 * t3426 - t3069 * t3422 - t3071 * t3421 - t3073 * t3420) * t3707 + (-t3098 * t3324 - t3111 * t3318 - t3112 * t3312 - t3113 * t3306) * MDP(12) + (t3098 * t3325 + t3111 * t3319 + t3112 * t3313 + t3113 * t3307) * MDP(13)) * t3251) * t3253; t3039 + (t3046 * t3497 + t3062 * t3489 + t3064 * t3487 + t3066 * t3485) * MDP(5) + (-t3046 * t3172 - t3062 * t3178 - t3064 * t3179 - t3066 * t3180) * MDP(6) + (t3042 * t3376 + t3056 * t3368 + t3057 * t3366 + t3058 * t3364) * t3559 + (-t3172 * t3354 - t3178 * t3351 - t3179 * t3348 - t3180 * t3345) * t3558 + ((t3047 * t3612 + t3063 * t3607 + t3065 * t3605 + t3067 * t3603) * MDP(4) + (t3047 * t3498 + t3063 * t3490 + t3065 * t3488 + t3067 * t3486) * MDP(7) + (t3047 * t3377 + t3063 * t3369 + t3065 * t3367 + t3067 * t3365) * t3707 + ((-t3047 * t3651 - t3063 * t3644 - t3065 * t3643 - t3067 * t3642) * MDP(4) + (-t3098 * t3415 - t3111 * t3412 - t3112 * t3409 - t3113 * t3406) * MDP(5) + (t3098 * t3633 + t3111 * t3630 + t3112 * t3627 + t3113 * t3624) * MDP(6) + (-t3047 * t3531 - t3063 * t3527 - t3065 * t3526 - t3067 * t3525) * MDP(7) + (-t3047 * t3426 - t3063 * t3422 - t3065 * t3421 - t3067 * t3420) * t3707 + (-t3098 * t3326 - t3111 * t3320 - t3112 * t3314 - t3113 * t3308) * MDP(12) + (t3098 * t3327 + t3111 * t3321 + t3112 * t3315 + t3113 * t3309) * MDP(13)) * t3251) * t3253; t3041 + (t3096 * t3612 + t3103 * t3607 + t3105 * t3605 + t3107 * t3603 + (-t3096 * t3651 - t3103 * t3644 - t3105 * t3643 - t3107 * t3642) * t3251) * t3691 + (t3095 * t3497 + t3102 * t3489 + t3104 * t3487 + t3106 * t3485 + (-t3098 * t3322 - t3111 * t3316 - t3112 * t3310 - t3113 * t3304) * t3565) * MDP(5) + (-t3095 * t3172 - t3102 * t3178 - t3104 * t3179 - t3106 * t3180 + t3268 * t3565) * MDP(6) + (t3096 * t3498 + t3103 * t3490 + t3105 * t3488 + t3107 * t3486 + (-t3096 * t3531 - t3103 * t3527 - t3105 * t3526 - t3107 * t3525) * t3251) * t3690 + (t3096 * t3377 + t3103 * t3369 + t3105 * t3367 + t3107 * t3365 + (-t3096 * t3426 - t3103 * t3422 - t3105 * t3421 - t3107 * t3420) * t3251) * t3557 + (-t3268 * t3252 + (t3200 * t3476 + t3207 * t3454 + t3211 * t3450 + t3215 * t3446) * t3251) * t3688 + ((-t3642 - t3643 - t3644 - t3651) * t3252 + (t3603 + t3605 + t3607 + t3612) * t3251) * t3686 + (-t3078 * t3612 - t3079 * t3607 - t3080 * t3605 - t3081 * t3603 + (-t3098 * t3413 - t3111 * t3410 - t3112 * t3407 - t3113 * t3404) * t3565) * MDP(12) + (t3172 * t3556 + t3178 * t3555 + t3179 * t3554 + t3180 * t3553 + (t3098 * t3323 + t3111 * t3317 + t3112 * t3311 + t3113 * t3305) * t3565) * MDP(13); (t3074 * t3612 + t3075 * t3607 + t3076 * t3605 + t3077 * t3603) * t3693 + (t3025 * t3612 + t3030 * t3607 + t3032 * t3605 + t3034 * t3603 + (-t3025 * t3651 - t3030 * t3644 - t3032 * t3643 - t3034 * t3642) * t3251) * t3691 + (t3024 * t3497 + t3029 * t3489 + t3031 * t3487 + t3033 * t3485 + (-t3098 * t3547 - t3111 * t3543 - t3112 * t3539 - t3113 * t3535) * t3251) * MDP(5) + (-t3024 * t3172 - t3029 * t3178 - t3031 * t3179 - t3033 * t3180 + (t3098 * t3548 + t3111 * t3544 + t3112 * t3540 + t3113 * t3536) * t3251) * MDP(6) + (t3025 * t3498 + t3030 * t3490 + t3032 * t3488 + t3034 * t3486 + (-t3025 * t3531 - t3030 * t3527 - t3032 * t3526 - t3034 * t3525) * t3251) * t3690 + (t3025 * t3377 + t3030 * t3369 + t3032 * t3367 + t3034 * t3365 + (-t3025 * t3426 - t3030 * t3422 - t3032 * t3421 - t3034 * t3420) * t3251) * t3557 + (t3172 * t3343 + t3178 * t3342 + t3179 * t3341 + t3180 * t3340 + (-t3098 * t3434 - t3111 * t3432 - t3112 * t3430 - t3113 * t3428) * t3251) * MDP(12) + (t3172 * t3335 + t3178 * t3334 + t3179 * t3333 + t3180 * t3332 + (t3098 * t3435 + t3111 * t3433 + t3112 * t3431 + t3113 * t3429) * t3251) * MDP(13) + t3291; t3039 + (t3048 * t3499 + t3068 * t3495 + t3070 * t3493 + t3072 * t3491) * MDP(5) + (-t3048 * t3171 - t3068 * t3175 - t3070 * t3176 - t3072 * t3177) * MDP(6) + (t3043 * t3378 + t3059 * t3374 + t3060 * t3372 + t3061 * t3370) * t3559 + (-t3171 * t3353 - t3175 * t3350 - t3176 * t3347 - t3177 * t3344) * t3558 + ((t3049 * t3613 + t3069 * t3610 + t3071 * t3609 + t3073 * t3608) * MDP(4) + (t3049 * t3500 + t3069 * t3496 + t3071 * t3494 + t3073 * t3492) * MDP(7) + (t3049 * t3379 + t3069 * t3375 + t3071 * t3373 + t3073 * t3371) * t3707 + ((t3049 * t3652 + t3069 * t3647 + t3071 * t3646 + t3073 * t3645) * MDP(4) + (t3097 * t3414 + t3108 * t3411 + t3109 * t3408 + t3110 * t3405) * MDP(5) + (-t3097 * t3632 - t3108 * t3629 - t3109 * t3626 - t3110 * t3623) * MDP(6) + (t3049 * t3532 + t3069 * t3530 + t3071 * t3529 + t3073 * t3528) * MDP(7) + (t3049 * t3427 + t3069 * t3425 + t3071 * t3424 + t3073 * t3423) * t3707 + (t3097 * t3324 + t3108 * t3318 + t3109 * t3312 + t3110 * t3306) * MDP(12) + (-t3097 * t3325 - t3108 * t3319 - t3109 * t3313 - t3110 * t3307) * MDP(13)) * t3251) * t3253; (t3046 * t3499 + t3062 * t3495 + t3064 * t3493 + t3066 * t3491) * MDP(5) + (-t3046 * t3171 - t3062 * t3175 - t3064 * t3176 - t3066 * t3177) * MDP(6) + (t3042 * t3378 + t3056 * t3374 + t3057 * t3372 + t3058 * t3370) * t3559 + (-t3171 * t3354 - t3175 * t3351 - t3176 * t3348 - t3177 * t3345) * t3558 + t3136 + (t3171 ^ 2 * t3201 + t3175 ^ 2 * t3208 + t3176 ^ 2 * t3212 + t3177 ^ 2 * t3216) * t3692 + ((t3047 * t3613 + t3063 * t3610 + t3065 * t3609 + t3067 * t3608) * MDP(4) + (t3047 * t3500 + t3063 * t3496 + t3065 * t3494 + t3067 * t3492) * MDP(7) + (t3047 * t3379 + t3063 * t3375 + t3065 * t3373 + t3067 * t3371) * t3707 + ((t3047 * t3652 + t3063 * t3647 + t3065 * t3646 + t3067 * t3645) * MDP(4) + (t3097 * t3415 + t3108 * t3412 + t3109 * t3409 + t3110 * t3406) * MDP(5) + (-t3097 * t3633 - t3108 * t3630 - t3109 * t3627 - t3110 * t3624) * MDP(6) + (t3047 * t3532 + t3063 * t3530 + t3065 * t3529 + t3067 * t3528) * MDP(7) + (t3047 * t3427 + t3063 * t3425 + t3065 * t3424 + t3067 * t3423) * t3707 + (t3097 * t3326 + t3108 * t3320 + t3109 * t3314 + t3110 * t3308) * MDP(12) + (-t3097 * t3327 - t3108 * t3321 - t3109 * t3315 - t3110 * t3309) * MDP(13)) * t3251) * t3253; t3040 + (t3096 * t3613 + t3103 * t3610 + t3105 * t3609 + t3107 * t3608 + (t3096 * t3652 + t3103 * t3647 + t3105 * t3646 + t3107 * t3645) * t3251) * t3691 + (t3095 * t3499 + t3102 * t3495 + t3104 * t3493 + t3106 * t3491 + (t3097 * t3322 + t3108 * t3316 + t3109 * t3310 + t3110 * t3304) * t3565) * MDP(5) + (-t3095 * t3171 - t3102 * t3175 - t3104 * t3176 - t3106 * t3177 - t3269 * t3565) * MDP(6) + (t3096 * t3500 + t3103 * t3496 + t3105 * t3494 + t3107 * t3492 + (t3096 * t3532 + t3103 * t3530 + t3105 * t3529 + t3107 * t3528) * t3251) * t3690 + (t3096 * t3379 + t3103 * t3375 + t3105 * t3373 + t3107 * t3371 + (t3096 * t3427 + t3103 * t3425 + t3105 * t3424 + t3107 * t3423) * t3251) * t3557 + (t3269 * t3252 + (t3200 * t3477 + t3207 * t3455 + t3211 * t3451 + t3215 * t3447) * t3251) * t3688 + ((t3645 + t3646 + t3647 + t3652) * t3252 + (t3608 + t3609 + t3610 + t3613) * t3251) * t3686 + (-t3078 * t3613 - t3079 * t3610 - t3080 * t3609 - t3081 * t3608 + (t3097 * t3413 + t3108 * t3410 + t3109 * t3407 + t3110 * t3404) * t3565) * MDP(12) + (t3171 * t3556 + t3175 * t3555 + t3176 * t3554 + t3177 * t3553 + (-t3097 * t3323 - t3108 * t3317 - t3109 * t3311 - t3110 * t3305) * t3565) * MDP(13); (t3074 * t3613 + t3075 * t3610 + t3076 * t3609 + t3077 * t3608) * t3693 + (t3025 * t3613 + t3030 * t3610 + t3032 * t3609 + t3034 * t3608 + (t3025 * t3652 + t3030 * t3647 + t3032 * t3646 + t3034 * t3645) * t3251) * t3691 + (t3024 * t3499 + t3029 * t3495 + t3031 * t3493 + t3033 * t3491 + (t3097 * t3547 + t3108 * t3543 + t3109 * t3539 + t3110 * t3535) * t3251) * MDP(5) + (-t3024 * t3171 - t3029 * t3175 - t3031 * t3176 - t3033 * t3177 + (-t3097 * t3548 - t3108 * t3544 - t3109 * t3540 - t3110 * t3536) * t3251) * MDP(6) + (t3025 * t3500 + t3030 * t3496 + t3032 * t3494 + t3034 * t3492 + (t3025 * t3532 + t3030 * t3530 + t3032 * t3529 + t3034 * t3528) * t3251) * t3690 + (t3025 * t3379 + t3030 * t3375 + t3032 * t3373 + t3034 * t3371 + (t3025 * t3427 + t3030 * t3425 + t3032 * t3424 + t3034 * t3423) * t3251) * t3557 + (t3171 * t3343 + t3175 * t3342 + t3176 * t3341 + t3177 * t3340 + (t3097 * t3434 + t3108 * t3432 + t3109 * t3430 + t3110 * t3428) * t3251) * MDP(12) + (t3171 * t3335 + t3175 * t3334 + t3176 * t3333 + t3177 * t3332 + (-t3097 * t3435 - t3108 * t3433 - t3109 * t3431 - t3110 * t3429) * t3251) * MDP(13) + t3440; t3041 + (t3049 * t3482 + t3069 * t3472 + t3071 * t3466 + t3073 * t3460 + (-t3049 * t3400 - t3069 * t3394 - t3071 * t3388 - t3073 * t3382) * t3251) * t3691 + (t3048 * t3362 + t3068 * t3360 + t3070 * t3358 + t3072 * t3356 + (-t3172 * t3301 - t3178 * t3298 - t3179 * t3295 - t3180 * t3292) * t3565) * MDP(5) + (-t3048 * t3597 - t3068 * t3583 - t3070 * t3581 - t3072 * t3579 + (t3172 * t3400 + t3178 * t3394 + t3179 * t3388 + t3180 * t3382) * t3565) * MDP(6) + (t3049 * t3484 + t3069 * t3474 + t3071 * t3468 + t3073 * t3462 + (-t3049 * t3403 - t3069 * t3397 - t3071 * t3391 - t3073 * t3385) * t3251) * t3690 + (t3049 * t3601 + t3069 * t3595 + t3071 * t3591 + t3073 * t3587 + (-t3049 * t3402 - t3069 * t3396 - t3071 * t3390 - t3073 * t3384) * t3251) * t3557 + (t3049 * t3597 + t3069 * t3583 + t3071 * t3581 + t3073 * t3579) * t3689 + (t3049 + t3069 + t3071 + t3073) * t3687 + (0.2e1 * t3353 + 0.2e1 * t3350 + 0.2e1 * t3347 + 0.2e1 * t3344 + (-t3476 - t3454 - t3450 - t3446 + (-t3172 * t3302 - t3178 * t3299 - t3179 * t3296 - t3180 * t3293) * t3253) * t3251) * MDP(12) + (t3043 * t3339 + t3059 * t3338 + t3060 * t3337 + t3061 * t3336 + (-t3172 - t3178 - t3179 - t3180 + (t3172 * t3303 + t3178 * t3300 + t3179 * t3297 + t3180 * t3294) * t3253) * t3251) * MDP(13); t3040 + (t3047 * t3482 + t3063 * t3472 + t3065 * t3466 + t3067 * t3460 + (-t3047 * t3400 - t3063 * t3394 - t3065 * t3388 - t3067 * t3382) * t3251) * t3691 + (t3046 * t3362 + t3062 * t3360 + t3064 * t3358 + t3066 * t3356 + (-t3171 * t3301 - t3175 * t3298 - t3176 * t3295 - t3177 * t3292) * t3565) * MDP(5) + (-t3046 * t3597 - t3062 * t3583 - t3064 * t3581 - t3066 * t3579 + (t3171 * t3400 + t3175 * t3394 + t3176 * t3388 + t3177 * t3382) * t3565) * MDP(6) + (t3047 * t3484 + t3063 * t3474 + t3065 * t3468 + t3067 * t3462 + (-t3047 * t3403 - t3063 * t3397 - t3065 * t3391 - t3067 * t3385) * t3251) * t3690 + (t3047 * t3601 + t3063 * t3595 + t3065 * t3591 + t3067 * t3587 + (-t3047 * t3402 - t3063 * t3396 - t3065 * t3390 - t3067 * t3384) * t3251) * t3557 + (t3047 * t3597 + t3063 * t3583 + t3065 * t3581 + t3067 * t3579) * t3689 + (t3047 + t3063 + t3065 + t3067) * t3687 + (0.2e1 * t3354 + 0.2e1 * t3351 + 0.2e1 * t3348 + 0.2e1 * t3345 + (-t3477 - t3455 - t3451 - t3447 + (-t3171 * t3302 - t3175 * t3299 - t3176 * t3296 - t3177 * t3293) * t3253) * t3251) * MDP(12) + (t3042 * t3339 + t3056 * t3338 + t3057 * t3337 + t3058 * t3336 + (-t3171 - t3175 - t3176 - t3177 + (t3171 * t3303 + t3175 * t3300 + t3176 * t3297 + t3177 * t3294) * t3253) * t3251) * MDP(13); (t3201 * t3596 + t3208 * t3582 + t3212 * t3580 + t3216 * t3578) * t3692 + (t3200 * t3478 + t3207 * t3456 + t3211 * t3452 + t3215 * t3448 + (-t3096 * t3400 - t3103 * t3394 - t3105 * t3388 - t3107 * t3382) * t3251) * t3691 + (t3095 * t3362 + t3102 * t3360 + t3104 * t3358 + t3106 * t3356 + (-t3199 * t3399 - t3206 * t3393 - t3210 * t3387 - t3214 * t3381) * t3565) * MDP(5) + (-t3095 * t3597 - t3102 * t3583 - t3104 * t3581 - t3106 * t3579 + t3267 * t3565) * MDP(6) + (t3096 * t3484 + t3103 * t3474 + t3105 * t3468 + t3107 * t3462 + (-t3096 * t3403 - t3103 * t3397 - t3105 * t3391 - t3107 * t3385) * t3251) * t3690 + (t3096 * t3601 + t3103 * t3595 + t3105 * t3591 + t3107 * t3587 + (-t3096 * t3402 - t3103 * t3396 - t3105 * t3390 - t3107 * t3384) * t3251) * t3557 + (-t3267 * t3564 + (t3478 + t3456 + t3452 + t3448 + (t3200 * t3596 + t3207 * t3582 + t3211 * t3580 + t3215 * t3578) * t3253) * t3251) * MDP(9) + ((-t3382 - t3388 - t3394 - t3400) * t3564 + (t3096 + t3103 + t3105 + t3107 + (t3460 + t3466 + t3472 + t3482) * t3253) * t3251) * MDP(10) + (t3203 + t3218 + t3221 + t3224) * MDP(11) * t3252 + (-t3078 * t3482 - t3079 * t3472 - t3080 * t3466 - t3081 * t3460 + (-t3596 - t3582 - t3580 - t3578 + (-t3294 - t3297 - t3300 - t3303) * t3253) * t3251) * MDP(12) + (t3482 * t3706 + t3472 * t3705 + t3466 * t3704 + t3460 * t3703 + (-t3597 - t3583 - t3581 - t3579 + (t3198 * t3399 + t3205 * t3393 + t3209 * t3387 + t3213 * t3381) * t3253) * t3251) * MDP(13) + MDP(17); (t3074 * t3482 + t3075 * t3472 + t3076 * t3466 + t3077 * t3460) * t3693 + (t3025 * t3482 + t3030 * t3472 + t3032 * t3466 + t3034 * t3460 + (-t3025 * t3400 - t3030 * t3394 - t3032 * t3388 - t3034 * t3382) * t3251) * t3691 + (t3024 * t3362 + t3029 * t3360 + t3031 * t3358 + t3033 * t3356 + (-t3203 * t3331 - t3218 * t3330 - t3221 * t3329 - t3224 * t3328) * t3251) * MDP(5) + (-t3024 * t3597 - t3029 * t3583 - t3031 * t3581 - t3033 * t3579 + (t3074 * t3507 + t3075 * t3505 + t3076 * t3503 + t3077 * t3501) * t3251) * MDP(6) + (t3025 * t3484 + t3030 * t3474 + t3032 * t3468 + t3034 * t3462 + (-t3025 * t3403 - t3030 * t3397 - t3032 * t3391 - t3034 * t3385) * t3251) * t3690 + (t3025 * t3601 + t3030 * t3595 + t3032 * t3591 + t3034 * t3587 + (-t3025 * t3402 - t3030 * t3396 - t3032 * t3390 - t3034 * t3384) * t3251) * t3557 + (t3025 * t3597 + t3030 * t3583 + t3032 * t3581 + t3034 * t3579) * t3689 + (t3025 + t3030 + t3032 + t3034) * t3687 + (0.2e1 * t3355 + 0.2e1 * t3352 + 0.2e1 * t3349 + 0.2e1 * t3346 + (-t3202 * t3331 - t3217 * t3330 - t3220 * t3329 - t3223 * t3328 + (-t3074 * t3475 - t3075 * t3453 - t3076 * t3449 - t3077 * t3445) * pkin(1)) * t3251) * MDP(12) + (t3023 * t3339 + t3026 * t3338 + t3027 * t3337 + t3028 * t3336 + (t3200 * t3074 * t3401 + t3207 * t3075 * t3395 + t3211 * t3076 * t3389 + t3215 * t3077 * t3383 + (-t3655 - t3657 - t3659 - t3661) * pkin(1)) * t3251) * MDP(13); (t3035 * t3049 + t3036 * t3069 + t3037 * t3071 + t3038 * t3073) * MDP(4) + (-t3275 - t3279 - t3283 - t3287) * MDP(6) + (t3049 * t3671 + t3069 * t3668 + t3071 * t3665 + t3073 * t3662) * MDP(7) + (t3049 * t3552 + t3069 * t3551 + t3071 * t3550 + t3073 * t3549) * t3707 + (MDP(5) * t3275 + t3712 * (t3061 * t3699 + t3663)) * t3585 + (MDP(5) * t3279 + t3711 * (t3060 * t3700 + t3666)) * t3589 + (MDP(5) * t3283 + t3710 * (t3059 * t3701 + t3669)) * t3593 + (MDP(5) * t3287 + t3709 * (t3043 * t3702 + t3672)) * t3599 + (t3114 * t3611 + t3115 * t3606 + t3116 * t3604 + t3117 * t3602) * t3692 + ((t3049 * t3641 + t3069 * t3639 + t3071 * t3637 + t3073 * t3635) * MDP(4) + (t3049 * t3524 + t3069 * t3523 + t3071 * t3522 + t3073 * t3521) * MDP(7) + (t3049 * t3419 + t3069 * t3418 + t3071 * t3417 + t3073 * t3416) * t3707) * t3253 + t3291; (t3035 * t3047 + t3036 * t3063 + t3037 * t3065 + t3038 * t3067) * MDP(4) + (-t3276 - t3280 - t3284 - t3288) * MDP(6) + (t3047 * t3671 + t3063 * t3668 + t3065 * t3665 + t3067 * t3662) * MDP(7) + (t3047 * t3552 + t3063 * t3551 + t3065 * t3550 + t3067 * t3549) * t3707 + (MDP(5) * t3276 + t3712 * (t3058 * t3699 + t3664)) * t3585 + (MDP(5) * t3280 + t3711 * (t3057 * t3700 + t3667)) * t3589 + (MDP(5) * t3284 + t3710 * (t3056 * t3701 + t3670)) * t3593 + (MDP(5) * t3288 + t3709 * (t3042 * t3702 + t3673)) * t3599 + (t3171 * t3640 + t3175 * t3638 + t3176 * t3636 + t3177 * t3634) * t3692 + ((t3047 * t3641 + t3063 * t3639 + t3065 * t3637 + t3067 * t3635) * MDP(4) + (t3047 * t3524 + t3063 * t3523 + t3065 * t3522 + t3067 * t3521) * MDP(7) + (t3047 * t3419 + t3063 * t3418 + t3065 * t3417 + t3067 * t3416) * t3707) * t3253 + t3440; (t3035 * t3096 + t3036 * t3103 + t3037 * t3105 + t3038 * t3107) * MDP(4) + (-t3270 - t3648 - t3649 - t3650 - t3653) * MDP(6) + (t3096 * t3671 + t3103 * t3668 + t3105 * t3665 + t3107 * t3662) * MDP(7) + (t3096 * t3552 + t3103 * t3551 + t3105 * t3550 + t3107 * t3549) * t3707 + ((-MDP(12) * t3081 + MDP(13) * t3703) * t3117 + (MDP(5) * t3648 + (t3234 * MDP(12) + (-MDP(13) * t3214 + MDP(5) * t3234) * t3223) * t3038) * t3241) * t3215 + ((-MDP(12) * t3080 + MDP(13) * t3704) * t3116 + (MDP(5) * t3649 + (t3232 * MDP(12) + (-MDP(13) * t3210 + MDP(5) * t3232) * t3220) * t3037) * t3239) * t3211 + ((-MDP(12) * t3079 + MDP(13) * t3705) * t3115 + (MDP(5) * t3650 + (t3230 * MDP(12) + (-MDP(13) * t3206 + MDP(5) * t3230) * t3217) * t3036) * t3237) * t3207 + ((-MDP(12) * t3078 + MDP(13) * t3706) * t3114 + (MDP(5) * t3653 + (t3226 * MDP(12) + (-MDP(13) * t3199 + MDP(5) * t3226) * t3202) * t3035) * t3229) * t3200 + (t3579 * t3634 + t3581 * t3636 + t3583 * t3638 + t3597 * t3640) * t3692 + ((t3096 * t3641 + t3103 * t3639 + t3105 * t3637 + t3107 * t3635) * MDP(4) + (t3096 * t3524 + t3103 * t3523 + t3105 * t3522 + t3107 * t3521) * MDP(7) + (t3096 * t3419 + t3103 * t3418 + t3105 * t3417 + t3107 * t3416) * t3707) * t3253 + (t3270 * MDP(9) + (t3035 + t3036 + t3037 + t3038) * MDP(10) + ((t3114 * t3482 + t3115 * t3472 + t3116 * t3466 + t3117 * t3460) * MDP(9) + (t3635 + t3637 + t3639 + t3641) * MDP(10)) * t3253) * t3251; (t3025 * t3035 + t3030 * t3036 + t3032 * t3037 + t3034 * t3038) * MDP(4) + (t3585 * t3674 + t3589 * t3675 + t3593 * t3676 + t3599 * t3677) * MDP(5) + (-t3674 - t3675 - t3676 - t3677) * MDP(6) + (t3025 * t3671 + t3030 * t3668 + t3032 * t3665 + t3034 * t3662) * MDP(7) + (t3025 * t3552 + t3030 * t3551 + t3032 * t3550 + t3034 * t3549) * t3707 + (t3114 * t3439 + t3115 * t3438 + t3116 * t3437 + t3117 * t3436) * t3559 + (-t3114 * t3355 - t3115 * t3352 - t3116 * t3349 - t3117 * t3346) * t3558 + MDP(14) + ((t3074 * t3641 + t3075 * t3639 + t3076 * t3637 + t3077 * t3635) * MDP(1) + (t3025 * t3641 + t3030 * t3639 + t3032 * t3637 + t3034 * t3635) * MDP(4) + (t3025 * t3524 + t3030 * t3523 + t3032 * t3522 + t3034 * t3521) * MDP(7) + (t3025 * t3419 + t3030 * t3418 + t3032 * t3417 + t3034 * t3416) * t3707) * t3253 + ((t3035 * t3660 + t3036 * t3658 + t3037 * t3656 + t3038 * t3654) * MDP(5) + (-t3035 * t3661 - t3036 * t3659 - t3037 * t3657 - t3038 * t3655) * MDP(6) + (t3035 * t3545 + t3036 * t3541 + t3037 * t3537 + t3038 * t3533) * MDP(12) + (-t3035 * t3546 - t3036 * t3542 - t3037 * t3538 - t3038 * t3534) * MDP(13)) * pkin(1);];
%% Postprocessing: Reshape Output
% From vec2mat_4_matlab.m
res = [t1(1), t1(2), t1(3), t1(4); t1(5), t1(6), t1(7), t1(8); t1(9), t1(10), t1(11), t1(12); t1(13), t1(14), t1(15), t1(16);];
MMX  = res;