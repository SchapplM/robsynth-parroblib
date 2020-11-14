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
% Datum: 2020-08-06 19:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR12V2G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2P2A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:26:30
% EndTime: 2020-08-06 19:26:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (210->74), mult. (237->133), div. (18->6), fcn. (162->18), ass. (0->70)
t50 = (pkin(2) + pkin(3));
t90 = -2 * t50;
t38 = sin(qJ(1,3));
t44 = cos(qJ(1,3));
t49 = pkin(5) - pkin(6);
t56 = t38 * pkin(1) - t49 * t44;
t37 = sin(qJ(2,3));
t68 = t37 * qJ(3,3);
t62 = 0.2e1 * t68;
t89 = (t38 * t62 + t56) * t50;
t40 = sin(qJ(1,2));
t46 = cos(qJ(1,2));
t55 = t40 * pkin(1) - t49 * t46;
t39 = sin(qJ(2,2));
t67 = t39 * qJ(3,2);
t61 = 0.2e1 * t67;
t88 = (t40 * t61 + t55) * t50;
t42 = sin(qJ(1,1));
t48 = cos(qJ(1,1));
t54 = t42 * pkin(1) - t49 * t48;
t41 = sin(qJ(2,1));
t66 = t41 * qJ(3,1);
t60 = 0.2e1 * t66;
t87 = t50 * (t42 * t60 + t54);
t86 = (t38 * qJ(3,3) + t56 * t37) * qJ(3,3);
t43 = cos(qJ(2,3));
t71 = t50 * t43;
t85 = 0.1e1 / (pkin(1) + t68 + t71) / qJ(3,3);
t45 = cos(qJ(2,2));
t70 = t50 * t45;
t84 = 0.1e1 / (pkin(1) + t67 + t70) / qJ(3,2);
t47 = cos(qJ(2,1));
t69 = t50 * t47;
t83 = 0.1e1 / (pkin(1) + t66 + t69) / qJ(3,1);
t82 = (t40 * qJ(3,2) + t55 * t39) * qJ(3,2);
t81 = (qJ(3,3) + t50) * (-qJ(3,3) + t50);
t80 = (qJ(3,2) + t50) * (-qJ(3,2) + t50);
t79 = (qJ(3,1) + t50) * (-qJ(3,1) + t50);
t78 = (t42 * qJ(3,1) + t54 * t41) * qJ(3,1);
t77 = t38 * t49;
t76 = t40 * t49;
t75 = t42 * t49;
t13 = pkin(1) * t37 + qJ(3,3);
t74 = t50 * t13;
t14 = pkin(1) * t39 + qJ(3,2);
t73 = t50 * t14;
t15 = pkin(1) * t41 + qJ(3,1);
t72 = t50 * t15;
t65 = qJ(3,1) * t90;
t64 = qJ(3,2) * t90;
t63 = qJ(3,3) * t90;
t59 = t38 * t81;
t58 = t40 * t80;
t57 = t42 * t79;
t36 = legFrame(1,2);
t35 = legFrame(2,2);
t34 = legFrame(3,2);
t33 = t47 ^ 2;
t32 = t45 ^ 2;
t31 = t43 ^ 2;
t21 = cos(t36);
t20 = cos(t35);
t19 = cos(t34);
t18 = sin(t36);
t17 = sin(t35);
t16 = sin(t34);
t9 = pkin(1) * qJ(3,1) - t41 * t79;
t8 = pkin(1) * qJ(3,2) - t39 * t80;
t7 = pkin(1) * qJ(3,3) - t37 * t81;
t1 = [((t18 * t65 + t21 * t57) * t33 + (-t18 * t9 + t21 * t87) * t47 + t21 * t78 + t18 * t72) * t83, ((-t18 * t57 + t21 * t65) * t33 + (-t18 * t87 - t9 * t21) * t47 - t18 * t78 + t21 * t72) * t83, (t48 * t33 * t79 + ((pkin(1) + t60) * t48 + t75) * t69 + (t15 * t48 + t41 * t75) * qJ(3,1)) * t83; ((t17 * t64 + t20 * t58) * t32 + (-t17 * t8 + t20 * t88) * t45 + t20 * t82 + t17 * t73) * t84, ((-t17 * t58 + t20 * t64) * t32 + (-t17 * t88 - t8 * t20) * t45 - t17 * t82 + t20 * t73) * t84, (t46 * t32 * t80 + ((pkin(1) + t61) * t46 + t76) * t70 + (t14 * t46 + t39 * t76) * qJ(3,2)) * t84; ((t16 * t63 + t19 * t59) * t31 + (-t16 * t7 + t19 * t89) * t43 + t19 * t86 + t16 * t74) * t85, ((-t16 * t59 + t19 * t63) * t31 + (-t16 * t89 - t7 * t19) * t43 - t16 * t86 + t19 * t74) * t85, (t44 * t31 * t81 + ((pkin(1) + t62) * t44 + t77) * t71 + (t13 * t44 + t37 * t77) * qJ(3,3)) * t85;];
Jinv  = t1;
