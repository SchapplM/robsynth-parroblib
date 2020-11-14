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
% Datum: 2020-08-07 03:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V2G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G3P3A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:29:06
% EndTime: 2020-08-07 03:29:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (234->100), mult. (504->209), div. (9->3), fcn. (453->26), ass. (0->87)
t23 = cos(pkin(4));
t95 = pkin(1) * t23;
t27 = sin(qJ(3,3));
t94 = pkin(2) * t27;
t30 = sin(qJ(3,2));
t93 = pkin(2) * t30;
t33 = sin(qJ(3,1));
t92 = pkin(2) * t33;
t37 = cos(qJ(2,3));
t91 = pkin(2) * t37;
t40 = cos(qJ(2,2));
t90 = pkin(2) * t40;
t43 = cos(qJ(2,1));
t89 = pkin(2) * t43;
t36 = cos(qJ(3,3));
t88 = pkin(3) * t36 ^ 2;
t39 = cos(qJ(3,2));
t87 = pkin(3) * t39 ^ 2;
t42 = cos(qJ(3,1));
t86 = pkin(3) * t42 ^ 2;
t22 = sin(pkin(4));
t85 = pkin(6) * t22;
t84 = t27 * pkin(3);
t83 = t30 * pkin(3);
t82 = t33 * pkin(3);
t28 = sin(qJ(2,3));
t81 = t22 * t28;
t29 = sin(qJ(1,3));
t80 = t22 * t29;
t31 = sin(qJ(2,2));
t79 = t22 * t31;
t32 = sin(qJ(1,2));
t78 = t22 * t32;
t34 = sin(qJ(2,1));
t77 = t22 * t34;
t35 = sin(qJ(1,1));
t76 = t22 * t35;
t75 = t23 * t29;
t74 = t23 * t32;
t73 = t23 * t35;
t38 = cos(qJ(1,3));
t72 = t23 * t38;
t41 = cos(qJ(1,2));
t71 = t23 * t41;
t44 = cos(qJ(1,1));
t70 = t23 * t44;
t69 = t27 * t22;
t45 = pkin(8) + pkin(7);
t68 = t29 * t45;
t67 = t30 * t22;
t66 = t32 * t45;
t65 = t33 * t22;
t64 = t35 * t45;
t63 = t38 * t45;
t62 = t41 * t45;
t61 = t44 * t45;
t60 = t45 * t28;
t59 = t45 * t31;
t58 = t45 * t34;
t57 = t38 * t69;
t56 = t41 * t67;
t55 = t44 * t65;
t54 = -pkin(2) * t28 + t45 * t37;
t53 = -pkin(2) * t31 + t45 * t40;
t52 = -pkin(2) * t34 + t45 * t43;
t51 = pkin(3) * t29 * t69 + (-pkin(2) * t75 + t63) * t28 + (t38 * pkin(2) + t23 * t68) * t37;
t50 = pkin(3) * t32 * t67 + (t41 * pkin(2) + t23 * t66) * t40 + (-pkin(2) * t74 + t62) * t31;
t49 = pkin(3) * t35 * t65 + (-pkin(2) * t73 + t61) * t34 + (t44 * pkin(2) + t23 * t64) * t43;
t48 = t54 * t22 - t23 * t84;
t47 = t53 * t22 - t23 * t83;
t46 = t52 * t22 - t23 * t82;
t26 = legFrame(1,2);
t25 = legFrame(2,2);
t24 = legFrame(3,2);
t18 = cos(t26);
t17 = cos(t25);
t16 = cos(t24);
t15 = sin(t26);
t14 = sin(t25);
t13 = sin(t24);
t6 = -t34 * t73 + t44 * t43;
t5 = -t31 * t74 + t41 * t40;
t4 = -t28 * t75 + t38 * t37;
t3 = 0.1e1 / (-(t34 * t95 + t43 * t85) * t86 + (((-pkin(6) + t82) * t89 - pkin(6) * t58 + pkin(1) * t82) * t22 + t52 * t95) * t42 + pkin(2) * (pkin(1) + t58 + t89) * t65);
t2 = 0.1e1 / (-(t31 * t95 + t40 * t85) * t87 + (((-pkin(6) + t83) * t90 - pkin(6) * t59 + pkin(1) * t83) * t22 + t53 * t95) * t39 + pkin(2) * (pkin(1) + t59 + t90) * t67);
t1 = 0.1e1 / (-(t28 * t95 + t37 * t85) * t88 + (((-pkin(6) + t84) * t91 - pkin(6) * t60 + pkin(1) * t84) * t22 + t54 * t95) * t36 + pkin(2) * (pkin(1) + t60 + t91) * t69);
t7 = [(-(t15 * t77 + t18 * t6) * t86 + (t46 * t15 - t49 * t18) * t42 - (t23 * t15 + t18 * t76) * t92) * t3, ((t15 * t6 - t18 * t77) * t86 + (t49 * t15 + t46 * t18) * t42 + (t15 * t76 - t23 * t18) * t92) * t3, ((t34 * t70 + t35 * t43) * t86 + (-pkin(3) * t55 + (t35 * pkin(2) - t23 * t61) * t43 + t34 * (pkin(2) * t70 + t64)) * t42 - pkin(2) * t55) * t3; (-(t14 * t79 + t17 * t5) * t87 + (t47 * t14 - t50 * t17) * t39 - (t23 * t14 + t17 * t78) * t93) * t2, ((t14 * t5 - t17 * t79) * t87 + (t50 * t14 + t47 * t17) * t39 + (t14 * t78 - t23 * t17) * t93) * t2, ((t31 * t71 + t32 * t40) * t87 + (-pkin(3) * t56 + (t32 * pkin(2) - t23 * t62) * t40 + t31 * (pkin(2) * t71 + t66)) * t39 - pkin(2) * t56) * t2; (-(t13 * t81 + t16 * t4) * t88 + (t48 * t13 - t51 * t16) * t36 - (t23 * t13 + t16 * t80) * t94) * t1, ((t13 * t4 - t16 * t81) * t88 + (t51 * t13 + t48 * t16) * t36 + (t13 * t80 - t23 * t16) * t94) * t1, ((t28 * t72 + t29 * t37) * t88 + (-pkin(3) * t57 + (t29 * pkin(2) - t23 * t63) * t37 + t28 * (pkin(2) * t72 + t68)) * t36 - pkin(2) * t57) * t1;];
Jinv  = t7;
