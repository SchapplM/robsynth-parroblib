% Calculate vector of inverse dynamics forces for parallel robot
% P3PPR2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m [3x1]
%   mass of all robot links (including platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:31
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX = P3PPR2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:31:33
% EndTime: 2018-12-20 17:31:33
% DurationCPUTime: 0.21s
% Computational Cost: add. (243->94), mult. (427->131), div. (0->0), fcn. (290->8), ass. (0->56)
t51 = m(1) + m(2);
t61 = -m(2) + t51;
t60 = koppelP(1,1);
t59 = koppelP(2,1);
t58 = koppelP(3,1);
t57 = koppelP(1,2);
t56 = koppelP(2,2);
t55 = koppelP(3,2);
t54 = mrSges(3,1);
t53 = mrSges(3,2);
t52 = xP(3);
t50 = xDDP(1);
t49 = xDDP(2);
t48 = xDDP(3);
t47 = legFrame(1,3);
t46 = legFrame(2,3);
t45 = legFrame(3,3);
t44 = xDP(3) ^ 2;
t43 = cos(t52);
t42 = sin(t52);
t41 = cos(t47);
t40 = cos(t46);
t39 = cos(t45);
t38 = sin(t47);
t37 = sin(t46);
t36 = sin(t45);
t35 = t41 ^ 2;
t34 = t40 ^ 2;
t33 = t39 ^ 2;
t32 = t38 ^ 2;
t31 = t37 ^ 2;
t30 = t36 ^ 2;
t29 = t41 * g(1) + t38 * g(2);
t28 = t40 * g(1) + t37 * g(2);
t27 = t39 * g(1) + t36 * g(2);
t26 = -t38 * g(1) + t41 * g(2);
t25 = -t37 * g(1) + t40 * g(2);
t24 = -t36 * g(1) + t39 * g(2);
t23 = -t42 * t57 + t43 * t60;
t22 = -t42 * t56 + t43 * t59;
t21 = -t42 * t55 + t43 * t58;
t20 = t42 * t60 + t43 * t57;
t19 = t42 * t59 + t43 * t56;
t18 = t42 * t58 + t43 * t55;
t17 = -t42 * t53 + t43 * t54;
t16 = t42 * t54 + t43 * t53;
t15 = t61 * t41 * t38;
t14 = t61 * t40 * t37;
t13 = t61 * t39 * t36;
t12 = -t20 * t48 - t23 * t44 + t50;
t11 = -t19 * t48 - t22 * t44 + t50;
t10 = -t18 * t48 - t21 * t44 + t50;
t9 = -t44 * t20 + t23 * t48 + t49;
t8 = -t44 * t19 + t22 * t48 + t49;
t7 = -t44 * t18 + t21 * t48 + t49;
t1 = [t13 * t7 + t14 * t8 + t15 * t9 - t16 * t48 - t44 * t17 + (t50 - g(1)) * m(3) + (t33 * t10 + t34 * t11 + t35 * t12 - t39 * t27 - t40 * t28 - t41 * t29) * t51 + (t30 * t10 + t31 * t11 + t32 * t12 + t36 * t24 + t37 * t25 + t38 * t26) * m(2); t13 * t10 + t14 * t11 + t15 * t12 - t44 * t16 + t17 * t48 + (t49 - g(2)) * m(3) + (-t36 * t27 - t37 * t28 - t38 * t29 + t30 * t7 + t31 * t8 + t32 * t9) * t51 + (-t39 * t24 - t40 * t25 - t41 * t26 + t33 * t7 + t34 * t8 + t35 * t9) * m(2); -t16 * t50 + t17 * t49 + Ifges(3,3) * t48 - (-g(1) * t54 - g(2) * t53) * t42 + t43 * (g(1) * t53 - g(2) * t54) + ((t12 * t41 + t38 * t9 - t29) * (-t20 * t41 + t23 * t38) + (t11 * t40 + t37 * t8 - t28) * (-t19 * t40 + t22 * t37) + (t10 * t39 + t36 * t7 - t27) * (-t18 * t39 + t21 * t36)) * t51 + ((t12 * t38 - t41 * t9 + t26) * (-t20 * t38 - t23 * t41) + (t11 * t37 - t40 * t8 + t25) * (-t19 * t37 - t22 * t40) + (t10 * t36 - t39 * t7 + t24) * (-t18 * t36 - t21 * t39)) * m(2);];
tauX  = t1;
