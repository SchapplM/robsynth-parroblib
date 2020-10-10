% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-28 23:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V3G2P4A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(1,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G2P4A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G2P4A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G2P4A1_Jinv: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G2P4A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G2P4A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-28 23:12:50
% EndTime: 2020-09-28 23:12:50
% DurationCPUTime: 0.37s
% Computational Cost: add. (174->96), mult. (462->186), div. (0->0), fcn. (516->42), ass. (0->130)
t58 = legFrame(6,2);
t147 = sin(t58);
t59 = legFrame(5,2);
t146 = sin(t59);
t60 = legFrame(4,2);
t145 = sin(t60);
t61 = legFrame(3,2);
t144 = sin(t61);
t62 = legFrame(2,2);
t143 = sin(t62);
t63 = legFrame(1,2);
t142 = sin(t63);
t141 = cos(qJ(2,1));
t140 = cos(qJ(2,2));
t139 = cos(qJ(2,3));
t75 = xP(4);
t45 = sin(t75);
t82 = koppelP(6,2);
t138 = t45 * t82;
t83 = koppelP(5,2);
t137 = t45 * t83;
t84 = koppelP(4,2);
t136 = t45 * t84;
t85 = koppelP(3,2);
t135 = t45 * t85;
t86 = koppelP(2,2);
t134 = t45 * t86;
t87 = koppelP(1,2);
t133 = t45 * t87;
t88 = koppelP(6,1);
t132 = t45 * t88;
t89 = koppelP(5,1);
t131 = t45 * t89;
t90 = koppelP(4,1);
t130 = t45 * t90;
t91 = koppelP(3,1);
t129 = t45 * t91;
t92 = koppelP(2,1);
t128 = t45 * t92;
t93 = koppelP(1,1);
t127 = t45 * t93;
t48 = cos(t75);
t126 = t48 * t82;
t125 = t48 * t83;
t124 = t48 * t84;
t123 = t48 * t85;
t122 = t48 * t86;
t121 = t48 * t87;
t120 = t48 * t88;
t119 = t48 * t89;
t118 = t48 * t90;
t117 = t48 * t91;
t116 = t48 * t92;
t115 = t48 * t93;
t49 = sin(qJ(2,6));
t114 = t49 * sin(qJ(1,6));
t51 = sin(qJ(2,5));
t113 = t51 * sin(qJ(1,5));
t53 = sin(qJ(2,4));
t112 = t53 * sin(qJ(1,4));
t111 = cos(qJ(1,6)) * t49;
t110 = cos(qJ(1,5)) * t51;
t109 = cos(qJ(1,4)) * t53;
t64 = sin(qJ(2,3));
t108 = t64 * sin(qJ(1,3));
t66 = sin(qJ(2,2));
t107 = t66 * sin(qJ(1,2));
t68 = sin(qJ(2,1));
t106 = t68 * sin(qJ(1,1));
t105 = cos(qJ(1,3)) * t64;
t104 = cos(qJ(1,2)) * t66;
t103 = cos(qJ(1,1)) * t68;
t74 = xP(5);
t47 = cos(t74);
t76 = koppelP(6,3);
t102 = t76 * t47;
t77 = koppelP(5,3);
t101 = t77 * t47;
t78 = koppelP(4,3);
t100 = t78 * t47;
t79 = koppelP(3,3);
t99 = t79 * t47;
t80 = koppelP(2,3);
t98 = t80 * t47;
t81 = koppelP(1,3);
t97 = t81 * t47;
t96 = cos(qJ(2,4));
t95 = cos(qJ(2,5));
t94 = cos(qJ(2,6));
t73 = xP(6);
t46 = cos(t73);
t44 = sin(t74);
t43 = sin(t73);
t42 = cos(t63);
t41 = cos(t62);
t40 = cos(t61);
t39 = cos(t60);
t38 = cos(t59);
t37 = cos(t58);
t30 = t42 * t106 - t142 * t141;
t29 = t41 * t107 - t143 * t140;
t28 = t40 * t108 - t144 * t139;
t27 = t142 * t106 + t42 * t141;
t26 = t143 * t107 + t41 * t140;
t25 = t144 * t108 + t40 * t139;
t24 = t39 * t112 - t145 * t96;
t23 = t38 * t113 - t146 * t95;
t22 = t37 * t114 - t147 * t94;
t21 = t145 * t112 + t39 * t96;
t20 = t146 * t113 + t38 * t95;
t19 = t147 * t114 + t37 * t94;
t18 = t44 * t81 + (-t43 * t87 + t46 * t93) * t47;
t17 = t44 * t80 + (-t43 * t86 + t46 * t92) * t47;
t16 = t44 * t79 + (-t43 * t85 + t46 * t91) * t47;
t15 = t44 * t78 + (-t43 * t84 + t46 * t90) * t47;
t14 = t44 * t77 + (-t43 * t83 + t46 * t89) * t47;
t13 = t44 * t76 + (-t43 * t82 + t46 * t88) * t47;
t12 = (t44 * t118 - t136) * t46 + (-t44 * t124 - t130) * t43 - t48 * t100;
t11 = (t44 * t119 - t137) * t46 + (-t44 * t125 - t131) * t43 - t48 * t101;
t10 = (t44 * t120 - t138) * t46 + (-t44 * t126 - t132) * t43 - t48 * t102;
t9 = (-t44 * t121 - t127) * t43 + (t44 * t115 - t133) * t46 - t48 * t97;
t8 = (-t44 * t133 + t115) * t43 + (t44 * t127 + t121) * t46 - t45 * t97;
t7 = (-t44 * t122 - t128) * t43 + (t44 * t116 - t134) * t46 - t48 * t98;
t6 = (-t44 * t123 - t129) * t43 + (t44 * t117 - t135) * t46 - t48 * t99;
t5 = -t45 * t98 + (t44 * t128 + t122) * t46 + t43 * (-t44 * t134 + t116);
t4 = -t45 * t99 + (t44 * t129 + t123) * t46 + t43 * (-t44 * t135 + t117);
t3 = -t45 * t100 + (t44 * t130 + t124) * t46 + t43 * (-t44 * t136 + t118);
t2 = -t45 * t101 + (t44 * t131 + t125) * t46 + t43 * (-t44 * t137 + t119);
t1 = -t45 * t102 + (t44 * t132 + t126) * t46 + t43 * (-t44 * t138 + t120);
t31 = [t30, -t27, t103, t8 * t103 - t27 * t9, -t18 * t103 - t30 * t9, -t27 * t18 - t30 * t8; t29, -t26, t104, t5 * t104 - t26 * t7, -t17 * t104 - t29 * t7, -t26 * t17 - t29 * t5; t28, -t25, t105, t4 * t105 - t25 * t6, -t16 * t105 - t28 * t6, -t25 * t16 - t28 * t4; t24, -t21, t109, t3 * t109 - t12 * t21, -t15 * t109 - t24 * t12, -t21 * t15 - t24 * t3; t23, -t20, t110, -t11 * t20 + t2 * t110, -t23 * t11 - t14 * t110, -t20 * t14 - t23 * t2; t22, -t19, t111, t1 * t111 - t10 * t19, -t22 * t10 - t13 * t111, -t22 * t1 - t19 * t13;];
Jinv  = t31;
