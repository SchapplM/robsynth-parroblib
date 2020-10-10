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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR12V1G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2P2A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:27:41
% EndTime: 2020-08-06 18:27:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (66->37), mult. (90->75), div. (15->6), fcn. (93->18), ass. (0->44)
t29 = cos(qJ(3,3));
t46 = pkin(3) * t29;
t30 = cos(qJ(1,3));
t45 = pkin(3) * t30;
t31 = cos(qJ(3,2));
t44 = pkin(3) * t31;
t32 = cos(qJ(1,2));
t43 = pkin(3) * t32;
t33 = cos(qJ(3,1));
t42 = pkin(3) * t33;
t34 = cos(qJ(1,1));
t41 = pkin(3) * t34;
t23 = sin(qJ(3,3));
t7 = t23 * pkin(3) + qJ(2,3);
t4 = 0.1e1 / t7;
t40 = 0.1e1 / t23 * t4;
t25 = sin(qJ(3,2));
t8 = t25 * pkin(3) + qJ(2,2);
t5 = 0.1e1 / t8;
t39 = 0.1e1 / t25 * t5;
t27 = sin(qJ(3,1));
t9 = t27 * pkin(3) + qJ(2,1);
t6 = 0.1e1 / t9;
t38 = 0.1e1 / t27 * t6;
t37 = t29 * qJ(2,3);
t36 = t31 * qJ(2,2);
t35 = t33 * qJ(2,1);
t28 = sin(qJ(1,1));
t26 = sin(qJ(1,2));
t24 = sin(qJ(1,3));
t22 = legFrame(1,2);
t21 = legFrame(2,2);
t20 = legFrame(3,2);
t16 = pkin(1) + pkin(5) + pkin(6);
t15 = cos(t22);
t14 = cos(t21);
t13 = cos(t20);
t12 = sin(t22);
t11 = sin(t21);
t10 = sin(t20);
t3 = qJ(2,1) * t34 - t16 * t28;
t2 = qJ(2,2) * t32 - t16 * t26;
t1 = qJ(2,3) * t30 - t16 * t24;
t17 = [((t12 * t42 - t3 * t15) * t27 + (t33 - 0.1e1) * (t33 + 0.1e1) * t15 * t41 + t12 * t35) * t38, ((t3 * t12 + t15 * t42) * t27 + (-t33 ^ 2 + 0.1e1) * t12 * t41 + t15 * t35) * t38, (t16 * t34 + t28 * t9) * t6; ((t11 * t44 - t2 * t14) * t25 + (t31 - 0.1e1) * (t31 + 0.1e1) * t14 * t43 + t11 * t36) * t39, ((t2 * t11 + t14 * t44) * t25 + (-t31 ^ 2 + 0.1e1) * t11 * t43 + t14 * t36) * t39, (t16 * t32 + t26 * t8) * t5; ((-t1 * t13 + t10 * t46) * t23 + (t29 - 0.1e1) * (t29 + 0.1e1) * t13 * t45 + t10 * t37) * t40, ((t1 * t10 + t13 * t46) * t23 + (-t29 ^ 2 + 0.1e1) * t10 * t45 + t13 * t37) * t40, (t16 * t30 + t24 * t7) * t4;];
Jinv  = t17;
