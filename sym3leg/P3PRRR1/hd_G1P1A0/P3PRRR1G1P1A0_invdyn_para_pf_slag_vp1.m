% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:45
% EndTime: 2020-03-09 21:14:47
% DurationCPUTime: 1.37s
% Computational Cost: add. (10395->198), mult. (6124->378), div. (1542->5), fcn. (6672->18), ass. (0->162)
t122 = pkin(7) + qJ(2,3);
t111 = qJ(3,3) + t122;
t101 = cos(t111);
t105 = sin(t122);
t108 = cos(t122);
t98 = sin(t111);
t200 = 0.1e1 / (-t101 * t105 + t108 * t98);
t123 = pkin(7) + qJ(2,2);
t112 = qJ(3,2) + t123;
t102 = cos(t112);
t106 = sin(t123);
t109 = cos(t123);
t99 = sin(t112);
t199 = 0.1e1 / (-t102 * t106 + t109 * t99);
t124 = pkin(7) + qJ(2,1);
t113 = qJ(3,1) + t124;
t100 = sin(t113);
t103 = cos(t113);
t107 = sin(t124);
t110 = cos(t124);
t198 = 0.1e1 / (t100 * t110 - t103 * t107);
t137 = 0.1e1 / pkin(3);
t139 = 0.1e1 / pkin(2);
t169 = t137 * t139;
t162 = t200 * t169;
t133 = xDP(1);
t126 = legFrame(3,3);
t114 = sin(t126);
t117 = cos(t126);
t71 = t101 * t117 - t114 * t98;
t52 = -pkin(2) * (t105 * t114 - t108 * t117) + t71 * pkin(3);
t174 = t52 * t133;
t132 = xDP(2);
t70 = t101 * t114 + t117 * t98;
t49 = pkin(2) * (t105 * t117 + t108 * t114) + t70 * pkin(3);
t177 = t49 * t132;
t142 = (t174 + t177) * t162;
t170 = t133 * t139;
t171 = t132 * t139;
t187 = t200 * t71;
t188 = t200 * t70;
t34 = t170 * t187 + t171 * t188;
t16 = -t142 + t34;
t197 = t16 * t142;
t161 = t199 * t169;
t127 = legFrame(2,3);
t115 = sin(t127);
t118 = cos(t127);
t73 = t102 * t118 - t115 * t99;
t53 = -pkin(2) * (t106 * t115 - t109 * t118) + t73 * pkin(3);
t173 = t53 * t133;
t72 = t102 * t115 + t118 * t99;
t50 = pkin(2) * (t106 * t118 + t109 * t115) + t72 * pkin(3);
t176 = t50 * t132;
t141 = (t173 + t176) * t161;
t185 = t199 * t73;
t186 = t199 * t72;
t35 = t170 * t185 + t171 * t186;
t17 = -t141 + t35;
t196 = t17 * t141;
t160 = t198 * t169;
t128 = legFrame(1,3);
t116 = sin(t128);
t119 = cos(t128);
t75 = -t100 * t116 + t103 * t119;
t54 = -pkin(2) * (t107 * t116 - t110 * t119) + t75 * pkin(3);
t172 = t54 * t133;
t74 = t100 * t119 + t103 * t116;
t51 = pkin(2) * (t107 * t119 + t110 * t116) + t74 * pkin(3);
t175 = t51 * t132;
t140 = (t172 + t175) * t160;
t183 = t198 * t75;
t184 = t198 * t74;
t36 = t170 * t183 + t171 * t184;
t18 = -t140 + t36;
t195 = t18 * t140;
t194 = t49 * t200;
t193 = t50 * t199;
t192 = t51 * t198;
t191 = t52 * t200;
t190 = t53 * t199;
t189 = t54 * t198;
t181 = t137 * t198;
t120 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t88 = m(3) * t120 + Icges(3,3);
t180 = t137 * t88;
t179 = t139 * t200;
t178 = t139 * t199;
t76 = -rSges(3,1) * t108 - t105 * rSges(3,2);
t77 = rSges(3,1) * t105 - rSges(3,2) * t108;
t148 = (-t101 * t76 + t77 * t98) * pkin(2);
t40 = Icges(3,3) + (t120 + t148) * m(3);
t168 = t137 * t40 * t200;
t78 = -rSges(3,1) * t109 - t106 * rSges(3,2);
t79 = rSges(3,1) * t106 - rSges(3,2) * t109;
t147 = (-t102 * t78 + t79 * t99) * pkin(2);
t41 = Icges(3,3) + (t120 + t147) * m(3);
t167 = t137 * t41 * t199;
t80 = -rSges(3,1) * t110 - t107 * rSges(3,2);
t81 = rSges(3,1) * t107 - rSges(3,2) * t110;
t146 = (t100 * t81 - t103 * t80) * pkin(2);
t42 = Icges(3,3) + (t120 + t146) * m(3);
t166 = t42 * t181;
t165 = t200 * t180;
t164 = t199 * t180;
t163 = t198 * t180;
t159 = 0.2e1 * pkin(2) * pkin(3);
t158 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t82 = -g(1) * t114 + g(2) * t117;
t85 = g(1) * t117 + g(2) * t114;
t157 = -t101 * (rSges(3,1) * t82 - rSges(3,2) * t85) + t98 * (rSges(3,1) * t85 + rSges(3,2) * t82);
t156 = t101 * t77 + t76 * t98;
t83 = -g(1) * t115 + g(2) * t118;
t86 = g(1) * t118 + g(2) * t115;
t155 = -t102 * (rSges(3,1) * t83 - rSges(3,2) * t86) + t99 * (rSges(3,1) * t86 + rSges(3,2) * t83);
t154 = t102 * t79 + t78 * t99;
t84 = -g(1) * t116 + g(2) * t119;
t87 = g(1) * t119 + g(2) * t116;
t153 = t100 * (rSges(3,1) * t87 + rSges(3,2) * t84) - t103 * (rSges(3,1) * t84 - rSges(3,2) * t87);
t152 = t100 * t80 + t103 * t81;
t151 = t101 * t108 + t105 * t98;
t150 = t102 * t109 + t106 * t99;
t149 = t100 * t107 + t103 * t110;
t145 = t151 * pkin(2);
t144 = t150 * pkin(2);
t143 = t149 * pkin(2);
t138 = pkin(2) ^ 2;
t136 = pkin(3) ^ 2;
t131 = xDDP(1);
t130 = xDDP(2);
t104 = t138 + t120;
t39 = (t104 + 0.2e1 * t146) * m(3) + t158;
t38 = (t104 + 0.2e1 * t147) * m(3) + t158;
t37 = (t104 + 0.2e1 * t148) * m(3) + t158;
t30 = (-t54 * t163 + t42 * t183) * t139;
t29 = (-t51 * t163 + t42 * t184) * t139;
t28 = (-t53 * t164 + t41 * t185) * t139;
t27 = (-t50 * t164 + t41 * t186) * t139;
t26 = (-t52 * t165 + t40 * t187) * t139;
t25 = (-t49 * t165 + t40 * t188) * t139;
t24 = (-t54 * t166 + t39 * t183) * t139;
t23 = (-t51 * t166 + t39 * t184) * t139;
t22 = (-t53 * t167 + t38 * t185) * t139;
t21 = (-t50 * t167 + t38 * t186) * t139;
t20 = (-t52 * t168 + t37 * t187) * t139;
t19 = (-t49 * t168 + t37 * t188) * t139;
t15 = -(t172 / 0.2e1 + t175 / 0.2e1) * t160 + t36;
t14 = -(t173 / 0.2e1 + t176 / 0.2e1) * t161 + t35;
t13 = -(t174 / 0.2e1 + t177 / 0.2e1) * t162 + t34;
t12 = (-pkin(3) * t195 + (pkin(3) * t18 + t36 * t143) * t36) * t198 * t139;
t11 = (-pkin(3) * t196 + (pkin(3) * t17 + t35 * t144) * t35) * t178;
t10 = (-pkin(3) * t197 + (pkin(3) * t16 + t34 * t145) * t34) * t179;
t9 = (-(-t149 * t15 * t159 - t136 * t18 - t138 * t36) * t36 * t181 - (pkin(3) + t143) * t198 * t195) * t139;
t8 = ((-t150 * t14 * t159 - t136 * t17 - t138 * t35) * t137 * t35 + (pkin(3) + t144) * t196) * t178;
t7 = ((-t151 * t13 * t159 - t136 * t16 - t138 * t34) * t137 * t34 + (pkin(3) + t145) * t197) * t179;
t6 = t42 * t12 - t88 * t9 + (-t152 * pkin(2) * t36 ^ 2 + t153) * m(3);
t5 = t41 * t11 + t88 * t8 + (-t154 * pkin(2) * t35 ^ 2 + t155) * m(3);
t4 = t40 * t10 + t88 * t7 + (-t156 * pkin(2) * t34 ^ 2 + t157) * m(3);
t3 = t39 * t12 - t42 * t9 + ((-rSges(2,1) * t84 + rSges(2,2) * t87) * t110 + t107 * (rSges(2,1) * t87 + rSges(2,2) * t84)) * m(2) + ((-0.2e1 * t140 * t15 * t152 + t107 * t87 - t84 * t110) * pkin(2) + t153) * m(3);
t2 = t38 * t11 + t41 * t8 + ((-rSges(2,1) * t83 + rSges(2,2) * t86) * t109 + t106 * (rSges(2,1) * t86 + rSges(2,2) * t83)) * m(2) + ((-0.2e1 * t14 * t141 * t154 + t106 * t86 - t83 * t109) * pkin(2) + t155) * m(3);
t1 = t37 * t10 + t40 * t7 + ((-rSges(2,1) * t82 + rSges(2,2) * t85) * t108 + t105 * (rSges(2,1) * t85 + rSges(2,2) * t82)) * m(2) + ((-0.2e1 * t13 * t142 * t156 + t105 * t85 - t82 * t108) * pkin(2) + t157) * m(3);
t31 = [(-g(1) + t131) * m(4) + (t1 * t187 + t2 * t185 + t3 * t183 + (-t6 * t189 - t5 * t190 - t4 * t191) * t137 + (t20 * t187 + t22 * t185 + t24 * t183 + (-t30 * t189 - t28 * t190 - t26 * t191) * t137) * t131 + (t20 * t188 + t22 * t186 + t24 * t184 + (-t30 * t192 - t28 * t193 - t26 * t194) * t137) * t130) * t139; (-g(2) + t130) * m(4) + (t1 * t188 + t2 * t186 + t3 * t184 + (-t6 * t192 - t5 * t193 - t4 * t194) * t137 + (t19 * t187 + t21 * t185 + t23 * t183 + (-t29 * t189 - t27 * t190 - t25 * t191) * t137) * t131 + (t19 * t188 + t21 * t186 + t23 * t184 + (-t29 * t192 - t27 * t193 - t25 * t194) * t137) * t130) * t139; (-g(3) + xDDP(3)) * ((3 * m(1)) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4));];
tauX  = t31;
