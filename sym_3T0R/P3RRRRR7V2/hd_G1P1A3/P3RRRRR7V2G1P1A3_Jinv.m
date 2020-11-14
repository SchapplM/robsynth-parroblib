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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2020-08-07 09:41
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G1P1A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G1P1A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G1P1A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G1P1A3_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G1P1A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G1P1A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:40:49
% EndTime: 2020-08-07 09:40:50
% DurationCPUTime: 0.44s
% Computational Cost: add. (444->149), mult. (390->115), div. (27->8), fcn. (141->93), ass. (0->82)
t116 = 2 * pkin(1);
t90 = (pkin(3) ^ 2);
t92 = (pkin(2) ^ 2);
t114 = -2 * t90 - 2 * t92;
t113 = 0.1e1 / pkin(2) / pkin(3);
t64 = (legFrame(3,3) + qJ(1,3));
t46 = (qJ(3,3) + t64);
t29 = qJ(2,3) + t46;
t47 = (-qJ(3,3) + t64);
t30 = -qJ(2,3) + t47;
t112 = cos(t29) + cos(t30);
t65 = (legFrame(2,3) + qJ(1,2));
t52 = (qJ(3,2) + t65);
t35 = qJ(2,2) + t52;
t53 = (-qJ(3,2) + t65);
t36 = -qJ(2,2) + t53;
t111 = cos(t35) + cos(t36);
t66 = (legFrame(1,3) + qJ(1,1));
t58 = (qJ(3,1) + t66);
t41 = qJ(2,1) + t58;
t59 = (-qJ(3,1) + t66);
t42 = -qJ(2,1) + t59;
t110 = cos(t41) + cos(t42);
t49 = qJ(2,3) + t64;
t50 = -qJ(2,3) + t64;
t109 = cos(t49) + cos(t50);
t55 = qJ(2,2) + t65;
t56 = -qJ(2,2) + t65;
t108 = cos(t55) + cos(t56);
t61 = qJ(2,1) + t66;
t62 = -qJ(2,1) + t66;
t107 = cos(t61) + cos(t62);
t105 = 2 * pkin(3);
t71 = (pkin(5) + pkin(6) + pkin(7));
t104 = 2 * t71;
t88 = 2 * qJ(2,1);
t60 = t88 + t66;
t89 = -2 * qJ(2,1);
t63 = t89 + t66;
t85 = 2 * qJ(2,2);
t54 = t85 + t65;
t86 = -2 * qJ(2,2);
t57 = t86 + t65;
t82 = 2 * qJ(2,3);
t48 = t82 + t64;
t83 = -2 * qJ(2,3);
t51 = t83 + t64;
t103 = t113 / 0.2e1;
t102 = -0.2e1 * sin(t29) - 0.2e1 * sin(t30);
t101 = -0.2e1 * sin(t35) - 0.2e1 * sin(t36);
t100 = -0.2e1 * sin(t41) - 0.2e1 * sin(t42);
t67 = sin((qJ(2,3) + qJ(3,3)));
t75 = sin(qJ(3,3));
t76 = sin(qJ(2,3));
t81 = 2 * qJ(3,3);
t99 = 0.1e1 / (t75 * t116 + (sin((t81 + qJ(2,3))) - t76) * pkin(3) + (-sin((qJ(2,3) - qJ(3,3))) + t67) * pkin(2)) * t103;
t68 = sin((qJ(2,2) + qJ(3,2)));
t77 = sin(qJ(3,2));
t78 = sin(qJ(2,2));
t84 = 2 * qJ(3,2);
t98 = 0.1e1 / (t77 * t116 + (sin((t84 + qJ(2,2))) - t78) * pkin(3) + (-sin((qJ(2,2) - qJ(3,2))) + t68) * pkin(2)) * t103;
t69 = sin((qJ(2,1) + qJ(3,1)));
t79 = sin(qJ(3,1));
t80 = sin(qJ(2,1));
t87 = 2 * qJ(3,1);
t97 = 0.1e1 / (t79 * t116 + (sin((t87 + qJ(2,1))) - t80) * pkin(3) + (-sin((qJ(2,1) - qJ(3,1))) + t69) * pkin(2)) * t103;
t96 = -0.2e1 * sin(t49) - 0.2e1 * sin(t50);
t95 = -0.2e1 * sin(t55) - 0.2e1 * sin(t56);
t94 = -0.2e1 * sin(t61) - 0.2e1 * sin(t62);
t45 = t89 + t59;
t44 = t88 + t58;
t43 = -2 * qJ(3,1) + t63;
t40 = t87 + t60;
t39 = t86 + t53;
t38 = t85 + t52;
t37 = -2 * qJ(3,2) + t57;
t34 = t84 + t54;
t33 = t83 + t47;
t32 = t82 + t46;
t31 = -2 * qJ(3,3) + t51;
t28 = t81 + t48;
t1 = [(cos(t66) * t114 + (-cos(t63) - cos(t60)) * t92 + (-cos(t43) - cos(t40)) * t90 + (t71 * t100 - t110 * t116) * pkin(3) + (t71 * t94 - t107 * t116 + (-cos(t45) - cos(t44) - cos(t59) - cos(t58)) * t105) * pkin(2)) * t97, (sin(t66) * t114 + (-sin(t63) - sin(t60)) * t92 + (-sin(t43) - sin(t40)) * t90 + (pkin(1) * t100 + t110 * t104) * pkin(3) + (t107 * t104 + pkin(1) * t94 + (-sin(t45) - sin(t44) - sin(t59) - sin(t58)) * t105) * pkin(2)) * t97, (-t80 * pkin(2) - pkin(3) * t69) / t79 * t113; (cos(t65) * t114 + (-cos(t57) - cos(t54)) * t92 + (-cos(t37) - cos(t34)) * t90 + (t71 * t101 - t111 * t116) * pkin(3) + (t71 * t95 - t108 * t116 + (-cos(t39) - cos(t38) - cos(t53) - cos(t52)) * t105) * pkin(2)) * t98, (sin(t65) * t114 + (-sin(t57) - sin(t54)) * t92 + (-sin(t37) - sin(t34)) * t90 + (pkin(1) * t101 + t111 * t104) * pkin(3) + (t108 * t104 + pkin(1) * t95 + (-sin(t39) - sin(t38) - sin(t53) - sin(t52)) * t105) * pkin(2)) * t98, (-t78 * pkin(2) - pkin(3) * t68) / t77 * t113; (cos(t64) * t114 + (-cos(t51) - cos(t48)) * t92 + (-cos(t31) - cos(t28)) * t90 + (t71 * t102 - t112 * t116) * pkin(3) + (t71 * t96 - t109 * t116 + (-cos(t33) - cos(t32) - cos(t47) - cos(t46)) * t105) * pkin(2)) * t99, (sin(t64) * t114 + (-sin(t51) - sin(t48)) * t92 + (-sin(t31) - sin(t28)) * t90 + (pkin(1) * t102 + t112 * t104) * pkin(3) + (t109 * t104 + pkin(1) * t96 + (-sin(t33) - sin(t32) - sin(t47) - sin(t46)) * t105) * pkin(2)) * t99, (-t76 * pkin(2) - pkin(3) * t67) / t75 * t113;];
Jinv  = t1;
