% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 00:55
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V2G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G1P1A2_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 00:53:37
% EndTime: 2020-08-07 00:53:37
% DurationCPUTime: 0.48s
% Computational Cost: add. (366->131), mult. (648->214), div. (9->6), fcn. (567->26), ass. (0->110)
t67 = cos(pkin(4));
t128 = pkin(1) * t67;
t71 = sin(qJ(3,3));
t127 = pkin(3) * t71;
t74 = sin(qJ(3,2));
t126 = pkin(3) * t74;
t77 = sin(qJ(3,1));
t125 = pkin(3) * t77;
t80 = cos(qJ(3,3));
t124 = pkin(3) * t80;
t83 = cos(qJ(3,2));
t123 = pkin(3) * t83;
t86 = cos(qJ(3,1));
t122 = pkin(3) * t86;
t81 = cos(qJ(2,3));
t121 = pkin(6) * t81;
t84 = cos(qJ(2,2));
t120 = pkin(6) * t84;
t87 = cos(qJ(2,1));
t119 = pkin(6) * t87;
t60 = t80 ^ 2;
t118 = t60 * pkin(3);
t62 = t83 ^ 2;
t117 = t62 * pkin(3);
t64 = t86 ^ 2;
t116 = t64 * pkin(3);
t59 = t67 ^ 2;
t115 = (t59 - 0.1e1) * pkin(6);
t66 = sin(pkin(4));
t114 = t66 * t67;
t72 = sin(qJ(2,3));
t113 = t67 * t72;
t75 = sin(qJ(2,2));
t112 = t67 * t75;
t78 = sin(qJ(2,1));
t111 = t67 * t78;
t110 = t67 * t80;
t109 = t67 * t83;
t108 = t67 * t86;
t89 = pkin(8) + pkin(7);
t43 = t89 * t81;
t44 = t89 * t84;
t45 = t89 * t87;
t56 = t81 * pkin(2);
t28 = t89 * t72 + t56;
t57 = t84 * pkin(2);
t29 = t89 * t75 + t57;
t58 = t87 * pkin(2);
t30 = t89 * t78 + t58;
t25 = pkin(1) + t28;
t107 = pkin(2) * t71 * t25;
t26 = pkin(1) + t29;
t106 = pkin(2) * t74 * t26;
t27 = pkin(1) + t30;
t105 = pkin(2) * t77 * t27;
t104 = t66 * t113;
t103 = t66 * t112;
t102 = t66 * t111;
t101 = t71 * pkin(6) + pkin(3);
t100 = t74 * pkin(6) + pkin(3);
t99 = t77 * pkin(6) + pkin(3);
t37 = t56 + pkin(1);
t98 = -t28 * pkin(6) + t37 * t127;
t38 = t57 + pkin(1);
t97 = -t29 * pkin(6) + t38 * t126;
t39 = t58 + pkin(1);
t96 = -t30 * pkin(6) + t39 * t125;
t95 = -t71 * (t37 * t72 + t89 + (t72 * t124 - t43) * t81) * t66 + (t81 * t124 + t25) * t110;
t94 = -t74 * (t38 * t75 + t89 + (t75 * t123 - t44) * t84) * t66 + (t84 * t123 + t26) * t109;
t93 = -t77 * (t39 * t78 + t89 + (t78 * t122 - t45) * t87) * t66 + (t87 * t122 + t27) * t108;
t61 = t81 ^ 2;
t22 = t61 * pkin(2) + t72 * t43 - pkin(2);
t31 = -t72 * pkin(2) + t43;
t46 = t61 - 0.2e1;
t92 = ((t46 * t127 - pkin(6)) * t80 + t71 * t22) * t114 - t72 * (-t101 + t118) - (t31 * t80 + (t101 - 0.2e1 * t118) * t72) * t59;
t63 = t84 ^ 2;
t23 = t63 * pkin(2) + t75 * t44 - pkin(2);
t32 = -t75 * pkin(2) + t44;
t47 = t63 - 0.2e1;
t91 = ((t47 * t126 - pkin(6)) * t83 + t74 * t23) * t114 - t75 * (-t100 + t117) - (t32 * t83 + (t100 - 0.2e1 * t117) * t75) * t59;
t65 = t87 ^ 2;
t24 = t65 * pkin(2) + t78 * t45 - pkin(2);
t33 = -t78 * pkin(2) + t45;
t48 = t65 - 0.2e1;
t90 = ((t48 * t125 - pkin(6)) * t86 + t77 * t24) * t114 - t78 * (-t99 + t116) - (t33 * t86 + (t99 - 0.2e1 * t116) * t78) * t59;
t88 = cos(qJ(1,1));
t85 = cos(qJ(1,2));
t82 = cos(qJ(1,3));
t79 = sin(qJ(1,1));
t76 = sin(qJ(1,2));
t73 = sin(qJ(1,3));
t70 = legFrame(1,3);
t69 = legFrame(2,3);
t68 = legFrame(3,3);
t55 = cos(t70);
t54 = cos(t69);
t53 = cos(t68);
t52 = sin(t70);
t51 = sin(t69);
t50 = sin(t68);
t21 = -t52 * t79 + t88 * t55;
t20 = t88 * t52 + t79 * t55;
t19 = -t51 * t76 + t85 * t54;
t18 = t85 * t51 + t76 * t54;
t17 = -t50 * t73 + t82 * t53;
t16 = t82 * t50 + t73 * t53;
t3 = 0.1e1 / (pkin(1) * (t45 + (-pkin(2) - t122) * t78) * t108 + t66 * (-t116 * t119 + t96 * t86 + t105));
t2 = 0.1e1 / (pkin(1) * (t44 + (-pkin(2) - t123) * t75) * t109 + t66 * (-t117 * t120 + t97 * t83 + t106));
t1 = 0.1e1 / (pkin(1) * (t43 + (-pkin(2) - t124) * t72) * t110 + t66 * (-t118 * t121 + t98 * t80 + t107));
t4 = [(-t90 * t20 + t93 * t21) * t3, (t93 * t20 + t90 * t21) * t3, ((-t33 * t114 + t115) * t86 + (-pkin(6) * t102 - t24 * t59 + t27 * t87) * t77 + (-(t48 * t59 - t65 + 0.1e1) * t77 * t86 + (0.2e1 * t64 - 0.1e1) * t102) * pkin(3)) / (-(pkin(1) * t111 + t66 * t119) * t116 + (t33 * t128 + t66 * t96) * t86 + t66 * t105); (-t91 * t18 + t94 * t19) * t2, (t94 * t18 + t91 * t19) * t2, ((-t32 * t114 + t115) * t83 + (-pkin(6) * t103 - t23 * t59 + t26 * t84) * t74 + (-(t47 * t59 - t63 + 0.1e1) * t74 * t83 + (0.2e1 * t62 - 0.1e1) * t103) * pkin(3)) / (-(pkin(1) * t112 + t66 * t120) * t117 + (t32 * t128 + t66 * t97) * t83 + t66 * t106); (-t92 * t16 + t95 * t17) * t1, (t95 * t16 + t92 * t17) * t1, ((-t31 * t114 + t115) * t80 + (-pkin(6) * t104 - t22 * t59 + t25 * t81) * t71 + (-(t46 * t59 - t61 + 0.1e1) * t71 * t80 + (0.2e1 * t60 - 0.1e1) * t104) * pkin(3)) / (-(pkin(1) * t113 + t66 * t121) * t118 + (t31 * t128 + t66 * t98) * t80 + t66 * t107);];
Jinv  = t4;
