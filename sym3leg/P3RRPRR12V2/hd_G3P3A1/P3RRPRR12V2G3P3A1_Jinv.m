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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:33
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR12V2G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3P3A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:33:11
% EndTime: 2020-08-06 19:33:11
% DurationCPUTime: 0.25s
% Computational Cost: add. (210->74), mult. (237->133), div. (18->6), fcn. (162->18), ass. (0->73)
t53 = (pkin(2) + pkin(3));
t93 = -2 * t53;
t44 = sin(qJ(2,1));
t51 = cos(qJ(1,1));
t45 = sin(qJ(1,1));
t52 = pkin(5) - pkin(6);
t69 = pkin(1) * t51 + t45 * t52;
t92 = qJ(3,1) * (qJ(3,1) * t51 + t69 * t44);
t42 = sin(qJ(2,2));
t49 = cos(qJ(1,2));
t43 = sin(qJ(1,2));
t70 = pkin(1) * t49 + t43 * t52;
t91 = qJ(3,2) * (qJ(3,2) * t49 + t70 * t42);
t40 = sin(qJ(2,3));
t47 = cos(qJ(1,3));
t41 = sin(qJ(1,3));
t71 = pkin(1) * t47 + t41 * t52;
t90 = qJ(3,3) * (qJ(3,3) * t47 + t71 * t40);
t68 = t40 * qJ(3,3);
t46 = cos(qJ(2,3));
t74 = t53 * t46;
t89 = 0.1e1 / (pkin(1) + t68 + t74) / qJ(3,3);
t67 = t42 * qJ(3,2);
t48 = cos(qJ(2,2));
t73 = t53 * t48;
t88 = 0.1e1 / (pkin(1) + t67 + t73) / qJ(3,2);
t66 = t44 * qJ(3,1);
t50 = cos(qJ(2,1));
t72 = t53 * t50;
t87 = 0.1e1 / (pkin(1) + t66 + t72) / qJ(3,1);
t37 = legFrame(3,2);
t19 = sin(t37);
t86 = t19 * t53;
t38 = legFrame(2,2);
t20 = sin(t38);
t85 = t20 * t53;
t39 = legFrame(1,2);
t21 = sin(t39);
t84 = t21 * t53;
t22 = cos(t37);
t83 = t22 * t53;
t23 = cos(t38);
t82 = t23 * t53;
t24 = cos(t39);
t81 = t24 * t53;
t80 = (qJ(3,3) + t53) * (-qJ(3,3) + t53);
t79 = (qJ(3,2) + t53) * (-qJ(3,2) + t53);
t78 = (qJ(3,1) + t53) * (-qJ(3,1) + t53);
t77 = t52 * t47;
t76 = t52 * t49;
t75 = t52 * t51;
t65 = qJ(3,1) * t93;
t64 = qJ(3,2) * t93;
t63 = qJ(3,3) * t93;
t62 = 0.2e1 * t68;
t61 = 0.2e1 * t67;
t60 = 0.2e1 * t66;
t59 = t47 * t80;
t58 = t49 * t79;
t57 = t51 * t78;
t36 = t50 ^ 2;
t35 = t48 ^ 2;
t34 = t46 ^ 2;
t15 = pkin(1) * t44 + qJ(3,1);
t14 = pkin(1) * t42 + qJ(3,2);
t13 = pkin(1) * t40 + qJ(3,3);
t9 = pkin(1) * qJ(3,1) - t44 * t78;
t8 = pkin(1) * qJ(3,2) - t42 * t79;
t7 = pkin(1) * qJ(3,3) - t40 * t80;
t6 = t51 * t60 + t69;
t5 = t49 * t61 + t70;
t4 = t47 * t62 + t71;
t1 = [((t21 * t65 + t24 * t57) * t36 + (-t21 * t9 + t6 * t81) * t50 + t24 * t92 + t15 * t84) * t87, ((-t21 * t57 + t24 * t65) * t36 + (-t9 * t24 - t6 * t84) * t50 - t21 * t92 + t15 * t81) * t87, (-t45 * t36 * t78 - ((pkin(1) + t60) * t45 - t75) * t72 - (t15 * t45 - t44 * t75) * qJ(3,1)) * t87; ((t20 * t64 + t23 * t58) * t35 + (-t20 * t8 + t5 * t82) * t48 + t23 * t91 + t14 * t85) * t88, ((-t20 * t58 + t23 * t64) * t35 + (-t8 * t23 - t5 * t85) * t48 - t20 * t91 + t14 * t82) * t88, (-t43 * t35 * t79 - ((pkin(1) + t61) * t43 - t76) * t73 - (t14 * t43 - t42 * t76) * qJ(3,2)) * t88; ((t19 * t63 + t22 * t59) * t34 + (-t19 * t7 + t4 * t83) * t46 + t22 * t90 + t13 * t86) * t89, ((-t19 * t59 + t22 * t63) * t34 + (-t7 * t22 - t4 * t86) * t46 - t19 * t90 + t13 * t83) * t89, (-t41 * t34 * t80 - ((pkin(1) + t62) * t41 - t77) * t74 - (t13 * t41 - t40 * t77) * qJ(3,3)) * t89;];
Jinv  = t1;
