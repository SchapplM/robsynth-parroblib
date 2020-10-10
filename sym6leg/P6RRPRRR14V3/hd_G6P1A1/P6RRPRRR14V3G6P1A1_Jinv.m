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
% Datum: 2020-09-28 23:18
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V3G6P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(1,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G6P1A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G6P1A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G6P1A1_Jinv: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G6P1A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G6P1A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-28 23:18:21
% EndTime: 2020-09-28 23:18:21
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
t74 = xP(5);
t47 = cos(t74);
t76 = koppelP(6,3);
t126 = t47 * t76;
t77 = koppelP(5,3);
t125 = t47 * t77;
t78 = koppelP(4,3);
t124 = t47 * t78;
t79 = koppelP(3,3);
t123 = t47 * t79;
t80 = koppelP(2,3);
t122 = t47 * t80;
t81 = koppelP(1,3);
t121 = t47 * t81;
t48 = cos(t75);
t120 = t48 * t82;
t119 = t48 * t83;
t118 = t48 * t84;
t117 = t48 * t85;
t116 = t48 * t86;
t115 = t48 * t87;
t114 = t48 * t88;
t113 = t48 * t89;
t112 = t48 * t90;
t111 = t48 * t91;
t110 = t48 * t92;
t109 = t48 * t93;
t49 = sin(qJ(2,6));
t108 = t49 * sin(qJ(1,6));
t51 = sin(qJ(2,5));
t107 = t51 * sin(qJ(1,5));
t53 = sin(qJ(2,4));
t106 = t53 * sin(qJ(1,4));
t105 = cos(qJ(1,6)) * t49;
t104 = cos(qJ(1,5)) * t51;
t103 = cos(qJ(1,4)) * t53;
t64 = sin(qJ(2,3));
t102 = t64 * sin(qJ(1,3));
t66 = sin(qJ(2,2));
t101 = t66 * sin(qJ(1,2));
t68 = sin(qJ(2,1));
t100 = t68 * sin(qJ(1,1));
t99 = cos(qJ(1,3)) * t64;
t98 = cos(qJ(1,2)) * t66;
t97 = cos(qJ(1,1)) * t68;
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
t30 = t42 * t100 - t142 * t141;
t29 = t41 * t101 - t143 * t140;
t28 = t40 * t102 - t144 * t139;
t27 = t142 * t100 + t42 * t141;
t26 = t143 * t101 + t41 * t140;
t25 = t144 * t102 + t40 * t139;
t24 = t39 * t106 - t145 * t96;
t23 = t38 * t107 - t146 * t95;
t22 = t37 * t108 - t147 * t94;
t21 = t145 * t106 + t39 * t96;
t20 = t146 * t107 + t38 * t95;
t19 = t147 * t108 + t37 * t94;
t18 = t44 * t81 + (-t43 * t87 + t46 * t93) * t47;
t17 = t44 * t80 + (-t43 * t86 + t46 * t92) * t47;
t16 = t44 * t79 + (-t43 * t85 + t46 * t91) * t47;
t15 = t44 * t78 + (-t43 * t84 + t46 * t90) * t47;
t14 = t44 * t77 + (-t43 * t83 + t46 * t89) * t47;
t13 = t44 * t76 + (-t43 * t82 + t46 * t88) * t47;
t12 = (t44 * t112 - t136) * t46 + (-t44 * t118 - t130) * t43 - t48 * t124;
t11 = (t44 * t113 - t137) * t46 + (-t44 * t119 - t131) * t43 - t48 * t125;
t10 = (t44 * t114 - t138) * t46 + (-t44 * t120 - t132) * t43 - t48 * t126;
t9 = (-t44 * t115 - t127) * t43 + (t44 * t109 - t133) * t46 - t48 * t121;
t8 = (-t44 * t133 + t109) * t43 + (t44 * t127 + t115) * t46 - t45 * t121;
t7 = (-t44 * t116 - t128) * t43 + (t44 * t110 - t134) * t46 - t48 * t122;
t6 = (-t44 * t117 - t129) * t43 + (t44 * t111 - t135) * t46 - t48 * t123;
t5 = -t45 * t122 + (t44 * t128 + t116) * t46 + t43 * (-t44 * t134 + t110);
t4 = -t45 * t123 + (t44 * t129 + t117) * t46 + t43 * (-t44 * t135 + t111);
t3 = -t45 * t124 + (t44 * t130 + t118) * t46 + t43 * (-t44 * t136 + t112);
t2 = -t45 * t125 + (t44 * t131 + t119) * t46 + t43 * (-t44 * t137 + t113);
t1 = -t45 * t126 + (t44 * t132 + t120) * t46 + t43 * (-t44 * t138 + t114);
t31 = [t30, -t27, t97, -t9 * t27 + t8 * t97, -t18 * t97 - t9 * t30, -t27 * t18 - t30 * t8; t29, -t26, t98, -t7 * t26 + t5 * t98, -t17 * t98 - t7 * t29, -t26 * t17 - t29 * t5; t28, -t25, t99, -t6 * t25 + t4 * t99, -t16 * t99 - t6 * t28, -t25 * t16 - t28 * t4; t24, -t21, t103, t3 * t103 - t12 * t21, -t15 * t103 - t12 * t24, -t21 * t15 - t24 * t3; t23, -t20, t104, t2 * t104 - t11 * t20, -t14 * t104 - t11 * t23, -t20 * t14 - t23 * t2; t22, -t19, t105, t1 * t105 - t10 * t19, -t10 * t22 - t13 * t105, -t22 * t1 - t19 * t13;];
Jinv  = t31;
